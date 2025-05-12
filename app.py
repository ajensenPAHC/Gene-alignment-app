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

            with open(fasta_file.name, 'r') as preview:
                fasta_preview = preview.read().split(">")
                if len(fasta_preview) > 1:
                    st.code(">" + fasta_preview[1], language="text")

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
