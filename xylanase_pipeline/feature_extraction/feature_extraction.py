# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 2025
@author: Bada Kanmi
"""
# xylanase_pipeline/feature_extraction/feature_extraction.py
from Bio import SeqIO  # For FASTA parsing
from Bio.SeqUtils.ProtParam import ProteinAnalysis  # For physicochemical props
import pandas as pd
import numpy as np
import os
from collections import Counter
import re  # Added for regex in extract_motifs

def load_fasta(fasta_path):
    """Load FASTA file into list of (header, sequence) tuples."""
    sequences = []
    with open(fasta_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            header = record.description  # Full header (e.g., ">P12345 | Name | Organism")
            seq = str(record.seq).upper()  # Clean uppercase sequence
            if len(seq) > 0:  # Skip empty
                sequences.append((header, seq))
    print(f"[INFO] Loaded {len(sequences)} sequences from {fasta_path}.")
    return sequences

def extract_aa_composition(seq):
    """Extract 20-dim AA composition (% frequency)."""
    aa_count = Counter(seq)
    total = len(seq)
    comp = {aa: (aa_count[aa] / total * 100) if total > 0 else 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}
    return pd.Series(comp)

def extract_physicochemical(seq):
    """Extract key physicochemical features via BioPython."""
    if len(seq) < 10:  # Skip too short
        return pd.Series({
            'length': len(seq),
            'molecular_weight': 0,
            'isoelectric_point': 0,
            'gravy': 0,  # Grand average of hydropathicity
            'instability_index': 0,
            'aromaticity': 0
        })
    analysis = ProteinAnalysis(seq)
    return pd.Series({
        'length': len(seq),
        'molecular_weight': analysis.molecular_weight(),
        'isoelectric_point': analysis.isoelectric_point(),
        'gravy': analysis.gravy(),  # Hydrophobicity
        'instability_index': analysis.instability_index(),
        'aromaticity': analysis.aromaticity()
    })

def extract_motifs(seq):
    """Basic motif detection for xylanases (e.g., catalytic residues)."""
    # Example: GH10/11 motifs (simplified regex for aspartate/glutamate pairs)
    # Basis: Conserved domains from CAZy/Pfam (e.g., GH10: [DE]x[DE] for active site)
    gh_motif = len(re.findall(r'[DE][A-Z]{0,2}[DE]', seq.upper()))  # Placeholder; customize
    return pd.Series({'gh_motif_count': gh_motif})

def extract_features(sequences):
    """Main: Extract all features into DataFrame."""
    features_list = []
    for header, seq in sequences:
        # Parse metadata from header (e.g., ">Acc | Name | Org")
        parts = header.split(' | ')
        acc = parts[0].replace('>', '') if len(parts) > 0 else 'Unknown'
        name = parts[1] if len(parts) > 1 else 'Unknown'
        org = ' | '.join(parts[2:]) if len(parts) > 2 else 'Unknown'
        
        # Extract features
        aa_comp = extract_aa_composition(seq)
        phys = extract_physicochemical(seq)
        motifs = extract_motifs(seq)
        
        row = pd.concat([aa_comp, phys, motifs])
        row['Accession'] = acc
        row['Protein_Name'] = name
        row['Organism'] = org
        row['Sequence'] = seq  # Optional: Keep raw seq
        features_list.append(row)
    
    df = pd.DataFrame(features_list)
    print(f"[INFO] Extracted {len(df)} feature vectors (shape: {df.shape}).")
    return df

def save_features(df, output_dir="../results/features", prefix="xylanase_features"):
    """Save features to CSV in output_dir."""
    os.makedirs(output_dir, exist_ok=True)
    csv_path = f"{output_dir}/{prefix}.csv"
    df.to_csv(csv_path, index=False)
    print(f"[INFO] Features saved to {csv_path}.")
    return csv_path

# Test block (conditional; skips if files missing)
if __name__ == "__main__":
    import os  # Ensure os is available
    # Test with fungal FASTA (absolute path; adjust if your username/path differs)
    fasta_path = "/home/kanmi77/projects/xylanase_pipeline/results/fasta/fungal_xylanase_sequences.fasta"
    if os.path.exists(fasta_path):
        sequences = load_fasta(fasta_path)
        df_features = extract_features(sequences)
        save_features(df_features, prefix="fungal_xylanase_features")
    else:
        print(f"[WARN] FASTA not found at {fasta_path}. Run retrieval pipeline first!")
        # Optional: Test with bacterial too
        bacterial_path = "/home/kanmi77/projects/xylanase_pipeline/results/fasta/bacterial_xylanase_sequences.fasta"
        if os.path.exists(bacterial_path):
            sequences = load_fasta(bacterial_path)
            df_features = extract_features(sequences)
            save_features(df_features, prefix="bacterial_xylanase_features")
        else:
            print("[INFO] No FASTA files found—extraction skipped. Run 'python -m sequence_retrieval.main' first.")
