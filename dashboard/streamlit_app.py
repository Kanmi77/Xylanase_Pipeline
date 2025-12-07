import streamlit as st
import pandas as pd
from pathlib import Path

st.set_page_config(page_title="Xylanase Dashboard", layout="wide")

ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = ROOT / 'results' / 'metadata'

st.title("Xylanase Sequences — Dashboard")
st.markdown("This dashboard loads metadata CSVs from the repository and provides interactive filtering and download.")

csv_files = list(METADATA_DIR.glob('*.csv'))
options = {f.name: f for f in csv_files}

choice = st.selectbox("Choose metadata file", options=list(options.keys()))

df = pd.read_csv(options[choice])

st.sidebar.header("Filters")
search = st.sidebar.text_input("Text search (any column)")

# Multiselect for Organism if present
if 'Organism' in df.columns:
    organisms = sorted(df['Organism'].dropna().unique())
    sel_org = st.sidebar.multiselect('Organism', organisms)
else:
    sel_org = []

if 'Sequence_Length' in df.columns:
    min_len = int(df['Sequence_Length'].min())
    max_len = int(df['Sequence_Length'].max())
    length_range = st.sidebar.slider('Sequence length', min_len, max_len, (min_len, max_len))
else:
    length_range = None

filtered = df.copy()
if search:
    mask = filtered.apply(lambda row: row.astype(str).str.contains(search, case=False).any(), axis=1)
    filtered = filtered[mask]
if sel_org:
    filtered = filtered[filtered['Organism'].isin(sel_org)]
if length_range is not None and 'Sequence_Length' in filtered.columns:
    filtered = filtered[filtered['Sequence_Length'].between(length_range[0], length_range[1])]

st.write(f"Showing {len(filtered)} rows — from `{choice}`")

st.dataframe(filtered)

csv_bytes = filtered.to_csv(index=False).encode('utf-8')
st.download_button('Download filtered CSV', data=csv_bytes, file_name='filtered_metadata.csv', mime='text/csv')

st.markdown("---")
st.markdown("Run locally: `streamlit run dashboard/streamlit_app.py`\n\nDeploy: use Streamlit Cloud or Render to host the app.")
