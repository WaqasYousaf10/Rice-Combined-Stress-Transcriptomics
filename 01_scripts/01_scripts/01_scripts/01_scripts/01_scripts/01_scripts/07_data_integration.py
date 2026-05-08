#!/usr/bin/env python3
# =============================================================================
# Data Integration: Physiology + Gene Expression
# =============================================================================
# Author: Waqas Yousaf
# Date: May 2026
# Usage: python 07_data_integration.py
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import os

# Create output directory
OUTPUT_DIR = "../results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 50)
print("DATA INTEGRATION PIPELINE")
print("=" * 50)

# Create example data (replace with actual data when available)
np.random.seed(42)
samples = []
genotypes = ["IR64", "IR64-Sub1"]
treatments = ["Control", "IronDeficiency", "Submergence", "Combined"]

for g in genotypes:
    for t in treatments:
        for r in [1, 2, 3]:
            samples.append(f"{g}_{t}_{r}")

# Create physiology data
physio_data = pd.DataFrame({
    "Sample_ID": samples,
    "Genotype": [s.split("_")[0] for s in samples],
    "Treatment": [s.split("_")[1] for s in samples],
    "Survival": np.random.uniform(20, 100, len(samples)),
    "Biomass": np.random.uniform(0.5, 5, len(samples)),
    "Fe_content": [150 if "Control" in s else 40 for s in samples],
    "ADH_activity": [10 if "Control" in s else 50 for s in samples],
})

# Create expression data
expr_data = pd.DataFrame({
    "Sample_ID": samples,
    "Genotype": [s.split("_")[0] for s in samples],
    "Treatment": [s.split("_")[1] for s in samples],
    "OsIRT1": np.random.uniform(0, 15, len(samples)),
    "SUB1A-1": [10 if "IR64-Sub1" in s else 1 for s in samples],
    "OsADH1": [5 if "Control" in s else 25 for s in samples],
})

# Merge datasets
merged_data = pd.merge(physio_data, expr_data, on=["Sample_ID", "Genotype", "Treatment"])

# Calculate correlations
correlations = []
for phys_var in ["Survival", "Biomass", "Fe_content"]:
    for expr_var in ["OsIRT1", "SUB1A-1", "OsADH1"]:
        corr, p_val = pearsonr(merged_data[phys_var], merged_data[expr_var])
        correlations.append({"Physiology": phys_var, "Gene": expr_var, 
                            "Correlation": corr, "P_value": p_val})

corr_df = pd.DataFrame(correlations)
corr_df.to_csv(os.path.join(OUTPUT_DIR, "correlations.csv"), index=False)

# Create correlation heatmap
pivot_corr = corr_df.pivot(index="Physiology", columns="Gene", values="Correlation")

plt.figure(figsize=(8, 6))
sns.heatmap(pivot_corr, annot=True, cmap="coolwarm", center=0, vmin=-1, vmax=1)
plt.title("Correlation: Physiology vs Gene Expression")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "correlation_heatmap.png"), dpi=300)
plt.close()

print(f"Results saved to: {OUTPUT_DIR}")
print("  - correlations.csv")
print("  - correlation_heatmap.png")
print("=" * 50)
