#!/bin/bash
#SBATCH --job-name=build_expression_matrix
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/build_expression_matrix.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/build_expression_matrix.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

# Load required modules
module load Python/3.12.3-GCCcore-13.3.0 Python-bundle-PyPI/2024.06-GCCcore-13.3.0

python3 - <<'EOF'
import pandas as pd
import numpy as np
from scipy.stats import norm as scipy_norm
from pathlib import Path

pd.set_option("future.no_silent_downcasting", True)

BASE = Path("/home/zw529/donglab/data/target_ALS")
OUT  = BASE / "eQTL"
OUT.mkdir(exist_ok=True)

def norm_id(x):
    return str(x).strip().replace("-", "_")

# 1. Load Metadata
rna = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")
rna.columns = rna.columns.str.strip().str.lower()
print(f"RNA samples in metadata: {len(rna)}")

# 2. Find normalization files (TPM tabs)
all_tabs = list(BASE.glob("**/RNAseq/Processed/*/normalization.tab"))
print(f"Found normalization.tab files: {len(all_tabs)}")

def find_tab(sample_id):
    sid_norm = norm_id(sample_id)
    for t in all_tabs:
        parent = norm_id(t.parent.name)
        if sid_norm == parent or sid_norm in parent:
            return t
    return None

rna["tab_path"] = rna["externalsampleid"].apply(find_tab)
rna_found = rna[rna["tab_path"].notna()].copy()
print(f"RNA samples with normalization.tab: {len(rna_found)}")

# 3. Build Expression Matrix (TPM)
expr = None
for _, row in rna_found.iterrows():
    sample = norm_id(row["externalsampleid"])
    tab = pd.read_csv(row["tab_path"], sep="\t")

    if not {"gene_id", "TPM"}.issubset(tab.columns):
        continue

    tab = tab[["gene_id", "TPM"]].copy()
    tab.columns = ["gene_id", sample]

    if expr is None:
        expr = tab
    else:
        expr = expr.merge(tab, on="gene_id", how="outer")

expr = expr.fillna(0).infer_objects(copy=False)
print(f"Raw matrix shape: {expr.shape}")

# 4. Rank-Based Inverse Normal Transformation (INT)
def inverse_normal_transform(series):
    n = len(series)
    # Rank data (averaging ties), map to (0, 1) range with offset
    ranks = series.rank(method='average')
    prob = (ranks - 0.5) / n
    # Map to standard normal Z-scores
    return scipy_norm.ppf(prob)

print("Applying Rank-Based Inverse Normal Transformation...")
gene_ids = expr['gene_id']
numeric_data = expr.drop('gene_id', axis=1)

# Apply INT per gene (row-wise)
# If a gene has 0 variance (all zeros), we keep it as 0 to avoid errors
transformed = numeric_data.apply(
    lambda x: inverse_normal_transform(x) if x.std() > 1e-9 else np.zeros(len(x)), 
    axis=1
)

# Reassemble
expr_final = pd.concat([gene_ids, transformed], axis=1)

# 5. Output
expr_final.to_csv(OUT / "expression_matrix.txt", sep="\t", index=False)

rna_found[["externalsampleid", "externalsubjectid", "tissue"]].to_csv(
    OUT / "expression_sample_metadata.csv", index=False
)

print(f"Saved normalized outputs to: {OUT}")
EOF
