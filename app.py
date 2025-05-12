import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import tempfile
import os
import subprocess
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

    row_mode = st.radio("How do you want to select rows?", ["Range", "Specific Indices"])
    if row_mode == "Range":
        row_range = st.slider("Select row range to process (based on index):", 0, len(df)-1, (1, len(df)-1))
        df = df.iloc[row_range[0]:row_range[1]+1]
    else:
        row_indices_input = st.text_input("Enter comma-separated row indices (e.g., 1,3,5):", "1,2")
        try:
            indices = [int(i.strip()) for i in row_indices_input.split(",") if i.strip().isdigit()]
            df = df.iloc[indices]
        except Exception as e:
            st.error(f"Invalid row indices: {e}")
            st.stop()

    sequences = []
    for index, row in df.iterrows():
        sequence = row[seq_column].strip().replace(" ", "") if pd.notna(row[seq_column]) else None
        if not sequence:
            continue
        name = "_".join([row[col].strip().replace(" ", "_") for col in name_columns if pd.notna(row[col])])
        if name and sequence:
            sequences.append(SeqRecord(Seq(sequence), id=name, description=""))

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
            fasta_path = fasta_file.name

        with open(fasta_path, 'r') as preview:
            fasta_preview = preview.read().split(">")
            if len(fasta_preview) > 1:
                st.code(">" + fasta_preview[1], language="text")

        with open(fasta_path, "rb") as download:
            st.download_button("Download FASTA Sequence File", download, file_name="sequences.fasta")

        # Run MAFFT locally
        aligned_path = fasta_path + ".aln"
        mafft_cline = MafftCommandline(input=fasta_path)
        stdout, stderr = mafft_cline()
        with open(aligned_path, "w") as aligned_file:
            aligned_file.write(stdout)

        alignment = AlignIO.read(aligned_path, "fasta")
        ref_aligned = next((rec for rec in alignment if rec.id == ref_seq.id), None)

        st.subheader("Visual Alignment Viewer")
        fig, ax = plt.subplots(figsize=(min(20, len(alignment[0].seq) / 4), len(alignment)))
        for i, record in enumerate(alignment):
            for j, (a, b) in enumerate(zip(record.seq, ref_aligned.seq)):
                match = a == b
                color = cm.Blues(1.0 if match else 0.0)
                ax.text(j, i, a, ha='center', va='center', fontsize=8, color='black',
                        bbox=dict(facecolor=color, edgecolor='none', boxstyle='round,pad=0.2'))
            ax.text(-1, i, record.id[:20], ha='right', va='center', fontsize=8)
        ax.set_xlim(-2, len(ref_aligned.seq))
        ax.set_ylim(-1, len(alignment))
        ax.axis('off')
        st.pyplot(fig)

        session["aligned_fasta"] = aligned_path
        session["ref_id"] = ref_seq.id
        st.success("MAFFT alignment complete!")
