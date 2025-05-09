import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests
import tempfile
import os
import time
from streamlit import components
from matplotlib import cm

st.set_page_config(page_title="Amino Acid Analyzer", layout="wide")
st.title("🧬 Amino Acid Sequence Analyzer and Classifier")

# Upload Excel
uploaded_excel = st.file_uploader("Upload Excel Spreadsheet (.xlsx)", type=["xlsx"])
uploaded_gene_db = st.file_uploader("Upload Gene Type Database (.csv)", type=["csv"])
ref_seq = None

if uploaded_excel:
    df = pd.read_excel(uploaded_excel)
    st.subheader("Preview of Uploaded Excel")
    st.dataframe(df.head())

    # Parse sequences from Excel
    sequences = []
    names = []
    for index, row in df.iterrows():
        case_id = str(row[0]).strip().replace(" ", "_") if pd.notna(row[0]) else None
        farm_name = str(row[2]).strip().replace(" ", "_") if pd.notna(row[2]) else None
        sequence = str(row[5]).strip().replace(" ", "") if pd.notna(row[5]) else None

        if case_id and sequence:
            name = f"{case_id}_{farm_name if farm_name else 'UnknownFarm'}"
            sequences.append(SeqRecord(Seq(sequence), id=name, description=""))
            names.append(name)

    if len(sequences) == 0:
        st.error("No valid sequences found in Excel. Please ensure at least Case ID and amino acid sequence are present.")
        st.stop()

    # Reference sequence selection
    ref_option = st.radio("Is your reference sequence in the uploaded Excel file?", ["Yes", "No"])
    if ref_option == "Yes":
        ref_index = st.selectbox("Select the row index of the reference sequence:", range(len(sequences)))
        ref_seq = sequences[ref_index]
    else:
        uploaded_fasta = st.file_uploader("Upload Reference Sequence (FASTA format)", type=["fasta"])
        if uploaded_fasta:
            ref_seq = list(SeqIO.parse(uploaded_fasta, "fasta"))[0]

    if not ref_seq:
        st.warning("Please provide a valid reference sequence to proceed.")
        st.stop()

    # Button to continue to alignment
    if st.button("Submit Sequences for Alignment"):
        # Write all sequences including reference to FASTA file
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

        # Wait for alignment to finish
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
        with tempfile.NamedTemporaryFile(delete=False, mode='w+') as aln_file:
            aln_file.write(aln)
            aln_file.flush()
            aligned_seqs = list(SeqIO.parse(aln_file.name, "fasta"))

        st.success("Alignment complete! Proceed to next step.")

        # Alignment Visualization
        st.subheader("🔍 Visual Alignment Viewer")
        ref_aligned = next((s for s in aligned_seqs if s.id == ref_seq.id), None)
        if ref_aligned:
            alignment_display = "<style>pre { font-size: 12px; font-family: monospace; }</style><pre>"
            for record in aligned_seqs:
                line = f"{record.id[:15]:<17}"
                for a, b in zip(record.seq, ref_aligned.seq):
                    color = "#c8e6c9" if a == b else "#ffcdd2"
                    line += f'<span style="background-color:{color};">{a}</span>'
                alignment_display += line + "\n"
            alignment_display += "</pre>"
            st.markdown(alignment_display, unsafe_allow_html=True)

        # Amino acid selection
        st.header("Select Amino Acid Positions")
        aa_pos_input = st.text_input("Enter comma-separated amino acid positions to extract:", "1,4,25")
        if st.button("Extract and Compare Positions"):
            positions = [int(p.strip()) for p in aa_pos_input.split(',') if p.strip().isdigit()]
            if not ref_aligned:
                st.error("Reference sequence not found in alignment.")
                st.stop()

            results = []
            for record in aligned_seqs:
                if record.id == ref_aligned.id:
                    continue
                aa_data = {}
                for pos in positions:
                    if pos - 1 < len(ref_aligned.seq):
                        aa_data[f"Pos_{pos}"] = record.seq[pos - 1]
                identity = sum(1 for a, b in zip(record.seq, ref_aligned.seq) if a == b) / len(ref_aligned.seq) * 100
                results.append({"Sequence Name": record.id, **aa_data, "Identity %": identity})

            result_df = pd.DataFrame(results)

            # Gene typing
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

            # Color-coding confidence
            def color_confidence(val):
                green = cm.Greens(val / 100)
                red = cm.Reds(1 - val / 100)
                hex_color = f"background-color: rgba({int(255*red[0])},{int(255*green[1])},0,0.5)"
                return hex_color

            st.subheader("Final Result Table")
            styled_df = result_df.style.background_gradient(cmap='RdYlGn', subset=["Identity %"])
            st.dataframe(styled_df, use_container_width=True)

            csv = result_df.to_csv(index=False).encode('utf-8')
            st.download_button("Download Results as CSV", csv, "results.csv", "text/csv")
