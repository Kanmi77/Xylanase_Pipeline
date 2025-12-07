# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 19:39:32 2025

@author: Bada Kanmi
"""

# xylanase_pipeline/sequence_retrieval/utils.py

import pandas as pd

def categorize_by_temperature(df):
    """Categorize enzymes by optimum temperature."""
    bins = [0, 45, 65, 120]
    labels = ["Mesophilic", "Moderately Thermophilic", "Thermophilic"]

    if "Optimum_Temperature" in df.columns:
        try:
            df["Optimum_Temperature"] = pd.to_numeric(df["Optimum_Temperature"], errors="coerce")
            df["Thermo_Class"] = pd.cut(df["Optimum_Temperature"], bins=bins, labels=labels)
        except Exception as e:
            print(f"[WARN] Could not categorize temperatures: {e}")
    else:
        df["Thermo_Class"] = "Unknown"
    return df
