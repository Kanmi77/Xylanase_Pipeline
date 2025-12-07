# Xylanase_Pipeline
Xylanase Thermostability Prediction Pipeline — a computational framework for retrieving, engineering, and analyzing fungal and bacterial xylanases using modern bioinformatics and structural modeling tools.  This repository contains the full workflow for the automated retrieval, processing of GH10 and GH11 xylanase sequences from fungi and bacteria
The pipeline integrates:

🔬 1. Sequence Retrieval (UniProt API)

Automated download of reviewed xylanase sequences

Cleaning, parsing, and standardization of metadata

FASTA generation and metadata export

Optional extraction of optimum temperature and pH when available

⚙️ 2. Feature Extraction (BioPython-based)

Amino acid composition

Physicochemical properties (MW, pI, GRAVY, instability, aromaticity)

GH10 / GH11 catalytic motif detection

Creation of high-quality ML-ready feature tables

🧠 3. Machine Learning Preparation (future module)

Feature engineering

Thermostability classification and regression models

Bacterial vs. fungal comparative analysis
(Module placeholder for future commits.)

🧬 4. Structural Modeling (future module)

AlphaFold2 modeling of selected sequences

Rosetta ΔΔG calculation for mutation assessment

Stability engineering workflow

📁 5. Organized Results Storage

FASTA files

Clean metadata

Extracted feature datasets

Purpose

This project supports a Master’s thesis focused on understanding and predicting the thermostability of xylanases, enzymes essential for:

biofuel production

food and feed processing

paper and pulp industries

green biotechnology

By combining bioinformatics, machine learning, and structural biology, the pipeline aims to identify motifs and features that distinguish thermophilic, mesophilic, and psychrophilic xylanases across bacteria and fungi.

Key Skills Demonstrated

Computational biology & enzyme informatics

Data engineering (API retrieval, cleaning, wrangling)

Feature engineering & ML prep

BioPython workflow design

AlphaFold & Rosetta structural modeling (in future modules)

Python package structuring and modular pipeline development

## Dashboard

This repository contains two dashboard options:

- A static interactive table suitable for GitHub Pages at `docs/index.html`. It loads `results/metadata/xylanase_sequences_metadata.csv` and displays it with DataTables.
- A Python Streamlit app at `dashboard/streamlit_app.py` for richer filtering and downloads.

How to enable GitHub Pages (serve `docs/`):

1. Go to your repository Settings → Pages.
2. Select branch `main` and folder `/docs` and save.
3. Your site will be available at `https://<username>.github.io/<repo>/` after a short build.

Run the Streamlit app locally:

```bash
pip install -r requirements.txt
streamlit run dashboard/streamlit_app.py
```

