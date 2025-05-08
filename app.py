import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import requests
import tempfile
import os

st.set_page_config(page_title="Amino Acid Analyzer", layout="wide")
st.title("🧬 Amino Acid Sequence Analyzer and Classifier")

# File upload
uploaded_excel = st.file_uploader("Upload Excel Spreadsheet (.xlsx)", type=["xlsx"])
uploaded_reference = None
uploaded_gene_db = st.file_uploader("Upload Gene Type Database (.csv)", type=["csv"])

ref_seq = None

if uploaded_excel:
    df = pd.read_excel(uploaded_excel)

    # Extract data from Excel
    st.subheader("Preview of Uploaded Excel")
    st.dataframe(df.head())

    sequences = []
    names = []

    for index, row in df.iterrows():
        case_id = str(row[0])
        farm_name = str(row[1])
        raw_sequence = str(row[5])
        seq_name = f"{case_id}_{farm_name}"
        try:
            seq_obj = Seq(raw_sequence)
            sequences.append((seq_name, seq_obj))
            names.append(seq_name)
        except Exception as e:
            st.warning(f"Error processing sequence in row {index + 2}: {e}")

    ref_option = st.radio("Is your reference sequence in the uploaded Excel file?", ["Yes", "No"])

    if ref_option == "Yes":
        ref_index = st.selectbox("Select the row index of the reference sequence:", range(len(sequences)))
        reference_name, reference_seq = sequences[ref_index]
    else:
        uploaded_reference = st.file_uploader("Upload Reference Sequence (FASTA format)", type=["fasta"])
        if uploaded_reference:
            ref_record = next(SeqIO.parse(uploaded_reference, "fasta"))
            reference_seq = ref_record.seq
            reference_name = ref_record.id

    if not reference_seq:
        st.warning("Please provide a valid reference sequence to proceed.")
        st.stop()

    # Write sequences to temp file
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as fasta_file:
        for name, seq in sequences:
            fasta_file.write(f">{name}\n{str(seq)}\n")
        input_fasta_path = fasta_file.name

    # Run Clustal Omega via web API
    st.info("Sending sequences to Clustal Omega for alignment...")
    with open(input_fasta_path, 'rb') as f:
        response = requests.post(
            'https://www.ebi.ac.uk/Tools/services/rest/clustalo/run',
            data={'email': 'your.email@example.com', 'stype': 'protein'},
            files={'sequence': f}
        )
    job_id = response.text.strip()

    import time
    time.sleep(5)

    while True:
        status = requests.get(f'https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}').text
        if status == 'FINISHED':
            break
        elif status == 'ERROR':
            st.error("Clustal Omega job failed.")
            st.stop()
        time.sleep(3)

    alignment = requests.get(f'https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-fasta').text
    with open(input_fasta_path, 'w') as temp_out:
        temp_out.write(alignment)
    aligned_sequences = list(SeqIO.parse(input_fasta_path, "fasta"))

    st.subheader("Confidence Comparison to Reference")
    confidence_data = []
    for seq_record in aligned_sequences:
        match_count = sum(1 for a, b in zip(seq_record.seq, reference_seq) if a == b)
        total = len(reference_seq)
        confidence = round((match_count / total) * 100, 2)
        confidence_data.append({
            "Sequence Name": seq_record.id,
            "Confidence (%)": confidence
        })

    result_df = pd.DataFrame(confidence_data)

    if uploaded_gene_db:
        type_df = pd.read_csv(uploaded_gene_db)
        type_match = []

        for record in aligned_sequences:
            assigned_type = "Unknown"
            for _, row in type_df.iterrows():
                gene_type = row['GeneType']
                positions = [int(p) for p in str(row['Positions']).split(',')]
                expected = str(row['AminoAcids']).split(',')

                if all(str(record.seq[pos - 1]) == expected[i] for i, pos in enumerate(positions)):
                    assigned_type = gene_type
                    break
            type_match.append(assigned_type)
        result_df["Predicted Type"] = type_match

    st.dataframe(result_df)

    csv = result_df.to_csv(index=False).encode('utf-8')
    st.download_button("Download Results as CSV", csv, "results.csv", "text/csv")
