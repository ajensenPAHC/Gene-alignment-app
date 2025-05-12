import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests
import tempfile
import os
import time
from matplotlib import cm
from matplotlib import pyplot as plt
import re

st.set_page_config(page_title="Amino Acid Analyzer", layout="wide")
st.title("🧬 Amino Acid Sequence Analyzer and Classifier")

uploaded_excel = st.file_uploader("Upload Excel Spreadsheet (.xlsx)", type=["xlsx"])
uploaded_gene_db = st.file_uploader("Upload Gene Type Database (.csv)", type=["csv"])
ref_seq = None
session = st.session_state

if uploaded_excel:
    xls = pd.ExcelFile(uploaded_excel)
    tab_name = st.selectbox("Select the Excel Sheet to Use:", xls.sheet_names)
    df = pd.read_excel(xls, sheet_name=tab_name)
    df = df.astype(str)
    st.subheader("Preview of Selected Sheet")
    st.dataframe(df.head())

    name_columns = st.multiselect("Select column(s) to use for naming sequences:", df.columns)
    seq_column = st.selectbox("Select the column containing amino acid sequences:", df.columns)

    row_range = st.slider("Select row range to process (based on index):", 0, len(df)-1, (1, len(df)-1))
    df = df.iloc[row_range[0]:row_range[1]+1]

    sequences = []
    names = []
    for index, row in df.iterrows():
        sequence = row[seq_column].strip().replace(" ", "") if pd.notna(row[seq_column]) else None
        if not sequence:
            continue
        name = "_".join([row[col].strip().replace(" ", "_") for col in name_columns if pd.notna(row[col])])
        if name and sequence:
            sequences.append(SeqRecord(Seq(sequence), id=name, description=""))
            names.append(name)

    if not sequences:
        st.error("No valid sequences found. Please ensure the name and sequence columns are selected correctly.")
        st.stop()

    ref_option = st.radio("How do you want to provide the reference sequence?", ["Select from Excel", "Upload FASTA File"])
    if ref_option == "Select from Excel":
        ref_index = st.selectbox("Select the row index of the reference sequence:", df.index)
        ref_row = df.loc[ref_index]
        ref_name = "_".join([ref_row[col].strip().replace(" ", "_") for col in name_columns if pd.notna(ref_row[col])])
        ref_seq_text = ref_row[seq_column].strip().replace(" ", "")
        ref_seq = SeqRecord(Seq(ref_seq_text), id=ref_name, description="")
    else:
        uploaded_fasta = st.file_uploader("Upload Reference Sequence (FASTA format)", type=["fasta"])
        if uploaded_fasta:
            ref_seq = list(SeqIO.parse(uploaded_fasta, "fasta"))[0]

    if not ref_seq:
        st.warning("Please provide a valid reference sequence to proceed.")
        st.stop()

    if st.button("Submit Sequences for Alignment"):
        all_seqs = [ref_seq] + [s for s in sequences if s.id != ref_seq.id]
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as fasta_file:
            SeqIO.write(all_seqs, fasta_file.name, "fasta")
            fasta_file.flush()
            with open(fasta_file.name, 'rb') as f:
                st.info("Submitting alignment job to Clustal Omega Web API...")
                response = requests.post(
                    'https://www.ebi.ac.uk/Tools/services/rest/clustalo/run',
                    data={'email': 'your.email@example.com', 'stype': 'protein'},
                    files={'sequence': f}
                )
                job_id = response.text.strip()

        with st.spinner("Running alignment..."):
            time.sleep(5)
            while True:
                status = requests.get(f'https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}').text
                if status == 'FINISHED':
                    break
                elif status == 'ERROR':
                    st.error("Clustal Omega job failed.")
                    st.stop()
                time.sleep(3)

        aln = requests.get(f'https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-fasta').text
        session["aligned_fasta"] = aln
        session["ref_id"] = ref_seq.id
        st.success("Alignment complete!")

if "aligned_fasta" in session and "ref_id" in session:
    with tempfile.NamedTemporaryFile(delete=False, mode='w+', suffix='.fasta') as aln_file:
        aln_file.write(session["aligned_fasta"])
        aln_file.flush()
        aligned_seqs = list(SeqIO.parse(aln_file.name, "fasta"))

    ref_aligned = next((s for s in aligned_seqs if s.id == session["ref_id"]), None)
    st.subheader("🔍 Alignment Snapshot (Image)")

    if ref_aligned:
        fig, ax = plt.subplots(figsize=(min(20, len(ref_aligned.seq) / 5), len(aligned_seqs)))
        for i, record in enumerate(aligned_seqs):
            for j, (a, b) in enumerate(zip(record.seq, ref_aligned.seq)):
                color = 'green' if a == b else 'red'
                ax.text(j, i, a, ha='center', va='center', fontsize=8, color=color)
            ax.text(-1, i, record.id[:20], ha='right', va='center', fontsize=8)
        ax.set_xlim(-2, len(ref_aligned.seq))
        ax.set_ylim(-1, len(aligned_seqs))
        ax.axis('off')

        image_path = os.path.join(tempfile.gettempdir(), 'alignment_snapshot.png')
        fig.savefig(image_path, bbox_inches='tight')
        st.pyplot(fig)
        with open(image_path, "rb") as file:
            st.download_button("Download Alignment Image", file, file_name="alignment_snapshot.png")

    st.header("Select Amino Acid Positions or Ranges")
    aa_pos_input = st.text_input("Enter comma-separated positions or ranges (e.g., 1,4,10-15,25):", "1,4,25")
    if st.button("Extract and Compare Positions"):
        positions = []
        for item in aa_pos_input.split(','):
            if '-' in item:
                start, end = map(int, item.strip().split('-'))
                positions.extend(range(start, end + 1))
            elif item.strip().isdigit():
                positions.append(int(item.strip()))

        results = []
        for record in aligned_seqs:
            if record.id == ref_aligned.id:
                continue
            aa_data = {}
            for pos in positions:
                if pos - 1 < len(ref_aligned.seq):
                    aa_data[f"Pos_{pos}"] = record.seq[pos - 1]
            identity = sum(1 for a, b in zip(record.seq, ref_aligned.seq) if a == b) / len(ref_aligned.seq) * 100
            results.append({"Sequence Name": record.id, **aa_data, "Identity": identity})

        result_df = pd.DataFrame(results)

        if uploaded_gene_db:
            gene_df = pd.read_csv(uploaded_gene_db)
            types = []
            for row in result_df.itertuples(index=False):
                best_type = "Unknown"
                best_score = 0
                for _, g in gene_df.iterrows():
                    match_count = 0
                    expected = {k: g[k] for k in g.index if k.startswith("Pos_")}
                    for pos_key, expected_aa in expected.items():
                        if hasattr(row, pos_key) and getattr(row, pos_key) == expected_aa:
                            match_count += 1
                    percent = (match_count / len(expected)) * 100 if expected else 0
                    if percent > best_score:
                        best_score = percent
                        best_type = f"{g['Type']} ({int(percent)}%)"
                types.append(best_type)
            result_df["Assigned Type"] = types

        st.subheader("Final Result Table")
        try:
            styled_df = result_df.style.background_gradient(cmap='RdYlGn', subset=["Identity"])
            st.dataframe(styled_df, use_container_width=True)
        except KeyError:
            st.dataframe(result_df, use_container_width=True)

        csv = result_df.to_csv(index=False).encode('utf-8')
        st.download_button("Download Results as CSV", csv, "results.csv", "text/csv")
