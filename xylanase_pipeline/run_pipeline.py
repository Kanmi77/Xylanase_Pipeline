# -*- coding: utf-8 -*-
"""
Xylanase Pipeline Runner: Retrieval + Feature Extraction
Created on Thu Nov 14 2025
@author: Bada Kanmi
"""
import sys
import os
import time  # For timing runs
import importlib.util  # For dynamic imports

# Add project root to path for local imports
project_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, project_root)

def run_retrieval():
    """Run the sequence retrieval pipeline."""
    print("\n[PIPELINE] Starting retrieval...")
    start_time = time.time()
    
    # Dynamic import for sequence_retrieval.main
    seq_retrieval_path = os.path.join(project_root, 'sequence_retrieval', 'main.py')
    if not os.path.exists(seq_retrieval_path):
        raise FileNotFoundError(f"sequence_retrieval/main.py not found at {seq_retrieval_path}. Check structure.")
    
    spec = importlib.util.spec_from_file_location("sequence_retrieval.main", seq_retrieval_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.main()  # Call the main function
    
    elapsed = time.time() - start_time
    print(f"[PIPELINE] Retrieval completed in {elapsed/60:.1f} minutes.")

def run_feature_extraction():
    """Run feature extraction on generated FASTAs."""
    print("\n[PIPELINE] Starting feature extraction...")
    start_time = time.time()
    
    # Dynamic import for feature_extraction
    feat_ext_path = os.path.join(project_root, 'feature_extraction', 'feature_extraction.py')
    if not os.path.exists(feat_ext_path):
        raise FileNotFoundError(f"feature_extraction/feature_extraction.py not found at {feat_ext_path}.")
    
    spec = importlib.util.spec_from_file_location("feature_extraction", feat_ext_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    load_fasta = module.load_fasta
    extract_features = module.extract_features
    save_features = module.save_features
    
    fasta_dir = os.path.join(project_root, "results", "fasta")
    features_dir = os.path.join(project_root, "results", "features")
    os.makedirs(features_dir, exist_ok=True)
    
    print(f"[DEBUG] Looking for FASTAs in {fasta_dir}")
    
    taxa = ["fungal", "bacterial"]
    for taxon in taxa:
        fasta_path = os.path.join(fasta_dir, f"{taxon}_xylanase_sequences.fasta")
        print(f"[DEBUG] Checking {fasta_path}")
        if os.path.exists(fasta_path):
            print(f"[PIPELINE] Extracting features for {taxon}...")
            sequences = load_fasta(fasta_path)
            df_features = extract_features(sequences)
            save_path = save_features(df_features, features_dir, prefix=f"{taxon}_xylanase_features")
            print(f"[PIPELINE] {taxon.capitalize()} features saved: {save_path}")
        else:
            print(f"[WARN] FASTA missing for {taxon}: {fasta_path}")
    
    elapsed = time.time() - start_time
    print(f"[PIPELINE] Feature extraction completed in {elapsed:.1f} seconds.")

if __name__ == "__main__":
    print(f"[PIPELINE START] Full Xylanase Pipeline — {time.strftime('%Y-%m-%d %H:%M:%S')}")
    run_retrieval()  # Step 1: Fetch sequences
    run_feature_extraction()  # Step 2: Extract features
    print(f"[PIPELINE DONE] All steps complete. Check 'results/' for outputs.")

