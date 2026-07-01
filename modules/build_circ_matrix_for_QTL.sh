#!/bin/bash
#SBATCH --job-name=build_circ_matrix
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/build_circ_matrix.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/build_circ_matrix.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=160G

set -euo pipefail

# 1. Clear existing Python paths to prevent loading 3.10 libraries
unset PYTHONPATH

# 2. Load the specific 2024a toolchain and Python 3.12
module purge
module load Python/3.12.3-GCCcore-13.3.0
module load Python-bundle-PyPI/2024.06-GCCcore-13.3.0

python3 - <<'EOF'
import pandas as pd
import numpy as np
from pathlib import Path

pd.set_option("future.no_silent_downcasting", True)

BASE = Path("/home/zw529/donglab/data/target_ALS")
OUT  = BASE / "QTL"
OUT.mkdir(exist_ok=True)

def norm_id(x):
    return str(x).strip().replace("-", "_")

# 1. Load Metadata
rna = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")
rna.columns = rna.columns.str.strip().str.lower()
print(f"RNA samples in metadata: {len(rna)}")

# 2. Find circularRNA files
all_circ = list(BASE.glob("**/RNAseq/Processed/*/circularRNA_known_circ_percentage.txt"))
print(f"Found circularRNA files: {len(all_circ)}")

def find_circ(sample_id):
    sid_norm = norm_id(sample_id)
    for c in all_circ:
        parent = norm_id(c.parent.name)
        if sid_norm == parent or sid_norm in parent:
            return c
    return None

rna["circ_path"] = rna["externalsampleid"].apply(find_circ)
rna_found = rna[rna["circ_path"].notna()].copy()
print(f"RNA samples with circularRNA file: {len(rna_found)}")

# 3. Build Circular RNA Percentage Matrix
expr = None
for _, row in rna_found.iterrows():
    sample = norm_id(row["externalsampleid"])
    
    # Files are space/whitespace-delimited based on head output
    circ = pd.read_csv(row["circ_path"], sep=r"\s+")
    
    if "circ_percent" not in circ.columns:
        continue

    # Create a unique, reproducible identifier for each circRNA locus
    circ["circ_id"] = (
        circ["chrom"].astype(str) + ":" + 
        circ["start"].astype(str) + "-" + 
        circ["end"].astype(str) + ":" + 
        circ["strand"].astype(str)
    )
    
    circ = circ[["circ_id", "circ_percent"]].copy()
    circ.columns = ["circ_id", sample]
    
    if expr is None:
        expr = circ
    else:
        expr = expr.merge(circ, on="circ_id", how="outer")

# Fill unmapped/missing features across samples with 0% circular percentage
expr = expr.fillna(0).infer_objects(copy=False)
print(f"Raw circRNA matrix shape: {expr.shape}")

# 4. Variance-Based Filtering
print("Filtering out zero-variance and low-variance circRNA features...")
circ_ids = expr['circ_id']
numeric_data = expr.drop('circ_id', axis=1)

# Calculate standard deviation across samples for each circRNA (row-wise)
row_sds = numeric_data.std(axis=1)

# Keep features with variations across your cohort (drops uniformly 0.0 or flat 1.0 variants)
min_sd_threshold = 0.005
variant_mask = row_sds > min_sd_threshold

filtered_numeric = numeric_data[variant_mask].copy()
filtered_circ_ids = circ_ids[variant_mask]

original_count = len(numeric_data)
filtered_count = len(filtered_numeric)
dropped_count = original_count - filtered_count

print(f"Original features: {original_count}")
print(f"Features remaining after SD > {min_sd_threshold} filter: {filtered_count}")
print(f"Dropped flat/invariant lines: {dropped_count} ({dropped_count / original_count * 100:.2f}%)")

# Recombine preserving the matrix format
expr_final = pd.concat([filtered_circ_ids.reset_index(drop=True), filtered_numeric.reset_index(drop=True)], axis=1)

# 5. Output
expr_final.to_csv(OUT / "circ_matrix.txt", sep="\t", index=False)

rna_found[["externalsampleid", "externalsubjectid", "tissue"]].to_csv(
    OUT / "circ_sample_metadata.csv", index=False
)

print(f"Saved filtered circRNA matrix to: {OUT}")
EOF
