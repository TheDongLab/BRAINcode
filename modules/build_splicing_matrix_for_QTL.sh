#!/bin/bash
#SBATCH --job-name=build_splicing_matrix
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/build_splicing_matrix.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/build_splicing_matrix.err
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

# 1. Metadata
rna = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")
rna.columns = rna.columns.str.strip().str.lower()
print(f"RNA samples in metadata: {len(rna)}")

# 2. Find PSI files
all_psi = list(BASE.glob("**/leafcutter/psi/*.leafcutter.PSI.tsv"))
print(f"Found PSI files: {len(all_psi)}")

def find_psi(sample_id):
    sid_norm = norm_id(sample_id)
    for p in all_psi:
        fname_sample = norm_id(p.stem.split(".leafcutter")[0])
        if sid_norm == fname_sample or sid_norm in fname_sample:
            return p
    return None

rna["psi_path"] = rna["externalsampleid"].apply(find_psi)
rna_found = rna[rna["psi_path"].notna()].copy()
print(f"RNA samples with PSI: {len(rna_found)}")

# 3. Build Splicing Matrix
expr = None
for _, row in rna_found.iterrows():
    sample = norm_id(row["externalsampleid"])
    psi = pd.read_csv(row["psi_path"], sep="\t")

    if "PSI" not in psi.columns:
        continue

    # unique junction ID
    psi["junction_id"] = (
        psi["chrom"] + ":" + psi["strand"] + ":" +
        psi["start"].astype(str) + "-" + psi["end"].astype(str)
    )

    psi = psi[["junction_id", "PSI"]].copy()
    psi.columns = ["junction_id", sample]

    if expr is None:
        expr = psi
    else:
        expr = expr.merge(psi, on="junction_id", how="outer")

# IMPORTANT:
# Missing junction in a sample's refined LeafCutter PSI file stays NA.
# Do NOT convert missing junctions to PSI=0.
expr = expr.infer_objects(copy=False)

print(f"Raw splicing matrix shape: {expr.shape}")
print(f"Total missing PSI values: {expr.isna().sum().sum()}")

# 4. Variance-Based Filtering
print("Filtering out low-information junctions...")

junction_ids = expr["junction_id"]
numeric_data = expr.drop("junction_id", axis=1)

# Calculate variability using only samples where PSI was actually observed
row_sds = numeric_data.std(axis=1, skipna=True)

# Require at least a few genuine PSI measurements
n_observed = numeric_data.notna().sum(axis=1)

min_sd_threshold = 0.005
min_observed_samples = 5

variant_mask = (
    (row_sds > min_sd_threshold) &
    (n_observed >= min_observed_samples)
)

filtered_numeric = numeric_data[variant_mask].copy()
filtered_junction_ids = junction_ids[variant_mask]

original_count = len(numeric_data)
filtered_count = len(filtered_numeric)
dropped_count = original_count - filtered_count

print(f"Original features: {original_count}")
print(
    f"Features remaining after SD > {min_sd_threshold} "
    f"and >= {min_observed_samples} observed samples: {filtered_count}"
)
print(
    f"Dropped low-information lines: {dropped_count} "
    f"({dropped_count / original_count * 100:.2f}%)"
)

# Recombine preserving matrix format
expr_final = pd.concat(
    [
        filtered_junction_ids.reset_index(drop=True),
        filtered_numeric.reset_index(drop=True)
    ],
    axis=1
)

# 5. Output
expr_final.to_csv(
    OUT / "splicing_matrix.txt",
    sep="\t",
    index=False,
    na_rep="NA"
)

rna_found[
    ["externalsampleid", "externalsubjectid", "tissue"]
].to_csv(
    OUT / "splicing_sample_metadata.csv",
    index=False
)

print(f"Saved filtered raw splicing matrix to: {OUT}")
EOF
