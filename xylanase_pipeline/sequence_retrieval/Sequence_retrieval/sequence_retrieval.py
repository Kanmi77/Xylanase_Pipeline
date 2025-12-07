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
FIELDS = [
    "accession", "id", "protein_name", "organism_name", "ec_number",
    "sequence", "length"
    # Fixed: entry_name -> id (valid field); others confirmed valid
]
def fetch_uniprot_sequences(query, size=200, retries=3):
    """
    Fetch GH10/GH11 xylanase sequences and metadata from UniProt.
    Includes fallback handling for API query errors and connectivity issues.
    """
    base_fields = ",".join(FIELDS)
    safe_query = urllib.parse.quote(query, safe="():")
    url = f"{UNIPROT_API}?query={safe_query}&format=tsv&fields={base_fields}&size={size}"

    print(f"[DEBUG] Requesting UniProt API: {url}")

    for attempt in range(1, retries + 1):
        try:
            response = requests.get(url, timeout=30)

            if response.status_code == 200:
                data = StringIO(response.text)
                df = pd.read_csv(data, sep="\t")
                print(f"[INFO] Retrieved {len(df)} entries from UniProt on attempt {attempt}.")
                return df

            elif response.status_code == 400:
                print(f"[WARN] Query failed (400). Attempt {attempt}/{retries}")
                print(f"[DEBUG] Response: {response.text[:300]}...")

                # Try fallback query (simpler search)
                fallback_query = 'xylanase AND reviewed:true AND taxonomy_id:4751'
                safe_fallback = urllib.parse.quote(fallback_query, safe="():")
                fallback_url = f"{UNIPROT_API}?query={safe_fallback}&format=tsv&fields={base_fields}&size={size}"
                print(f"[DEBUG] Retrying with fallback query: {fallback_query}")
                response = requests.get(fallback_url, timeout=30)

                if response.status_code == 200:
                    data = StringIO(response.text)
                    df = pd.read_csv(data, sep="\t")
                    print(f"[INFO] Retrieved {len(df)} entries using fallback query.")
                    return df

                print(f"[ERROR] Fallback query failed with status {response.status_code}.")

            else:
                print(f"[WARN] Unexpected status code {response.status_code}. Attempt {attempt}/{retries}")

        except requests.exceptions.RequestException as e:
            print(f"[ERROR] Connection issue: {e}. Attempt {attempt}/{retries}")
        
        time.sleep(2)  # Wait a bit before retrying

    raise Exception("[FAIL] UniProt retrieval failed after multiple attempts.")

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
                    # Broader regex for varied phrasing
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
        "Sequence": "Sequence"
        # Fixed: Keys match actual TSV headers (capitalized); no CC fields
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
    df = df.drop_duplicates(subset=["Accession"])
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
            f.write(f">{row['Accession']} | {row['Protein_Name']} | {row['Organism']}\n{row['Sequence']}\n")
    print(f"[INFO] FASTA saved to {fasta_path}")
    return csv_path, fasta_path