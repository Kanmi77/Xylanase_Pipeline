# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 19:39:16 2025
@author: Bada Kanmi
"""
# xylanase_pipeline/sequence_retrieval/sequence_retrieval.py
import requests
import pandas as pd
from io import StringIO
from datetime import datetime
import os
import urllib.parse
import re  # For parsing in fetch_entry_details
import time  # For rate limiting

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
FIELDS_METADATA = [
    "accession", "id", "protein_name", "organism_name", "ec", "length"
    # ec confirmed valid
]
FIELDS_WITH_SEQ = FIELDS_METADATA + ["sequence"]

def fetch_uniprot_sequences(query, size=200, include_sequences=True):
    """Fetch GH10/GH11 xylanase sequences and metadata from UniProt."""
    fallback_needed = include_sequences and size > 25
    if fallback_needed:
        print(f"[INFO] size={size} > 25 with sequences; fetching metadata first, then sequences individually.")
        # Fetch metadata with large size
        df = _fetch_core(query, size, FIELDS_METADATA)
        # Quick rename for loop access (TSV uses "Entry", not "accession")
        if "Entry" in df.columns:
            df = df.rename(columns={"Entry": "accession"})
        # Fetch sequences individually
        df["Sequence"] = None
        for idx, acc in enumerate(df['accession']):
            seq_url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
            seq_resp = requests.get(seq_url)
            if seq_resp.status_code == 200:
                # Parse FASTA: Join non-header lines
                lines = [line.strip() for line in seq_resp.text.split('\n') if line.strip() and not line.startswith('>')]
                df.loc[df['accession'] == acc, 'Sequence'] = ''.join(lines)
            else:
                print(f"[WARN] Failed sequence fetch for {acc}: {seq_resp.status_code}")
            time.sleep(0.3)  # Rate limit
            if (idx + 1) % 20 == 0:  # Less frequent logging for larger sets
                print(f"[INFO] Fetched sequences for {idx+1}/{len(df)} entries")
        return df
    else:
        fields = FIELDS_WITH_SEQ if include_sequences else FIELDS_METADATA
        max_size = min(size, 25 if include_sequences else 500)
        return _fetch_core(query, max_size, fields)

def _fetch_core(query, size, fields):
    """Internal: Core fetch logic."""
    query_encoded = urllib.parse.quote(query)
    url = f"{UNIPROT_API}?query={query_encoded}&format=tsv&fields={','.join(fields)}&size={size}"
    # print(f"[DEBUG] Requesting: {url}")  # Uncomment if needed
    response = requests.get(url)
    if response.status_code != 200:
        print(f"[ERROR] Failed. Response body: {response.text}")
        raise Exception(f"Failed to retrieve data: {response.status_code}")
    data = StringIO(response.text)
    df = pd.read_csv(data, sep="\t")
    print(f"[INFO] Retrieved {len(df)} entries from UniProt.")
    return df

def fetch_entry_details(accession):
    """Fetch and parse optimum temperature/pH from full entry JSON."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        comments = data.get('comments', [])
        for comment in comments:
            if comment.get('commentType') == 'BIOPHYSICOCHEMICAL PROPERTIES':
                texts = comment.get('texts', [])
                if texts:
                    text = texts[0].get('value', '')
                    temp_match = re.search(r'(?:optimum|optimal)\s*(?:temperature|temp)[:\s.]*(\d+)[°\s]?(?:C|°C)?', text, re.I)
                    ph_match = re.search(r'(?:optimum|optimal)\s*pH[:\s.]*([\d.]+)', text, re.I)
                    return temp_match.group(1) if temp_match else None, ph_match.group(1) if ph_match else None
    return None, None

def clean_metadata(df):
    """Clean and rename UniProt columns for clarity."""
    rename_map = {
        "Entry": "Accession",
        "Entry Name": "ID",
        "Protein names": "Protein_Name",
        "Organism": "Organism",
        "EC number": "EC",
        "Length": "Sequence_Length",
        "Sequence": "Sequence",
        "accession": "Accession"  # Handle temp rename from fallback
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
    df = df.drop_duplicates(subset=["Accession"])
    # Filter non-null sequences if present
    if "Sequence" in df.columns:
        df = df[df["Sequence"].notnull()]
    print(f"[INFO] Cleaned dataset: {len(df)} unique sequences retained.")
    return df

def save_outputs(df, prefix="xylanase_sequences"):
    """Save metadata and FASTA files in results/ directories."""
    os.makedirs("../results/fasta", exist_ok=True)
    os.makedirs("../results/metadata", exist_ok=True)
    csv_path = f"../results/metadata/{prefix}_metadata.csv"
    fasta_path = f"../results/fasta/{prefix}.fasta"
    df.to_csv(csv_path, index=False)
    print(f"[INFO] Metadata saved to {csv_path}")
    with open(fasta_path, "w") as f:
        for _, row in df.iterrows():
            seq = row.get('Sequence', '')
            header = f">{row['Accession']} | {row.get('Protein_Name', '')} | {row.get('Organism', '')}"
            f.write(f"{header}\n{seq}\n")
    print(f"[INFO] FASTA saved to {fasta_path}")
    return csv_path, fasta_path
