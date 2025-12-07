# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 19:34:44 2025
@author: Bada Kanmi
"""
# xylanase_pipeline/sequence_retrieval/main.py
from sequence_retrieval.sequence_retrieval import (
    fetch_uniprot_sequences, clean_metadata, save_outputs, fetch_entry_details
)
from sequence_retrieval.utils import categorize_by_temperature
from datetime import datetime
import pandas as pd
import time

def process_taxon(query, taxon_name, size=200):
    """Process a single taxon: Fetch, clean, parse optima, categorize, save."""
    print(f"\n[START] Processing {taxon_name} — {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    df = fetch_uniprot_sequences(query=query, size=size, include_sequences=True)
    if len(df) == 0:
        print(f"[WARN] No results for {taxon_name}; skipping.")
        return
    df_clean = clean_metadata(df)
    
    # Add columns for optima
    df_clean["Optimum_Temperature"] = None
    df_clean["Optimum_pH"] = None
    
    # Parse optima
    print(f"[INFO] Parsing optimum temperature/pH for {taxon_name}...")
    total = len(df_clean)
    if total > 0:
        for idx, acc in enumerate(df_clean['Accession']):
            temp, ph = fetch_entry_details(acc)
            df_clean.loc[df_clean['Accession'] == acc, 'Optimum_Temperature'] = temp
            df_clean.loc[df_clean['Accession'] == acc, 'Optimum_pH'] = ph
            if (idx + 1) % 10 == 0 or (idx + 1) == total:
                print(f"[INFO] Processed {idx + 1}/{total} entries")
            time.sleep(0.3)
    else:
        print(f"[WARN] No entries to parse for {taxon_name}.")
    
    df_temp = categorize_by_temperature(df_clean)
    # Save with taxon-specific prefix
    save_outputs(df_temp, prefix=f"{taxon_name.lower()}_xylanase_sequences")
    print(f"[DONE] {taxon_name} processing completed.\n")
def main():
    print(f"\n[OVERALL START] Combined Bacterial + Fungal Retrieval — {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Base query components (shared)
    base_query = '(GH10 OR GH11) AND (xylanase OR "beta-xylosidase") AND reviewed:true'
    
    # Fungal query
    fungal_query = f'{base_query} AND taxonomy_id:4751'
    fungal_df = process_taxon(fungal_query, "Fungal", size=200)
    
    # Bacterial query
    bacterial_query = f'{base_query} AND taxonomy_id:2'
    bacterial_df = process_taxon(bacterial_query, "Bacterial", size=200)
    
    # Optional: Merge summaries (uncomment if you want a combined CSV)
    # if fungal_df is not None and bacterial_df is not None:
    #     fungal_df["Taxon"] = "Fungal"
    #     bacterial_df["Taxon"] = "Bacterial"
    #     combined = pd.concat([fungal_df, bacterial_df], ignore_index=True)
    #     combined.to_csv("../results/metadata/combined_xylanase_sequences_metadata.csv", index=False)
    #     print("[INFO] Combined summary saved to ../results/metadata/combined_xylanase_sequences_metadata.csv")
    
    print(f"[OVERALL DONE] Pipeline completed. Check ../results/ for separate fungal/bacterial files.\n")

if __name__ == "__main__":
    main()
