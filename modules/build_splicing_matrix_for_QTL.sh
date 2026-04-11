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
# We use the bundle which includes compatible numpy, scipy, and pandas
module purge
module load Python/3.12.3-GCCcore-13.3.0
module load Python-bundle-PyPI/2024.06-GCCcore-13.3.0

python3 - <<'EOF'
import pandas as pd
import numpy as np
from scipy.stats import norm as scipy_norm
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

expr = expr.fillna(0).infer_objects(copy=False)
print(f"Raw splicing matrix shape: {expr.shape}")

# 4. Rank-Based Inverse Normal Transformation (INT)
def inverse_normal_transform(series):
    n = len(series)
    # Rank data, map to (0, 1) probability space
    ranks = series.rank(method='average')
    prob = (ranks - 0.5) / n
    # Map to standard normal distribution
    return scipy_norm.ppf(prob)

print("Applying Rank-Based Inverse Normal Transformation to PSI values...")
junction_ids = expr['junction_id']
numeric_data = expr.drop('junction_id', axis=1)
sample_names = numeric_data.columns

transformed = numeric_data.apply(
    lambda x: inverse_normal_transform(x) if x.std() > 1e-9 else np.zeros(len(x)), 
    axis=1, 
    result_type='expand'
)

transformed.columns = sample_names

expr_final = pd.concat([junction_ids, transformed], axis=1)

# 5. Output
expr_final.to_csv(OUT / "splicing_matrix.txt", sep="\t", index=False)

rna_found[["externalsampleid", "externalsubjectid", "tissue"]].to_csv(
    OUT / "splicing_sample_metadata.csv", index=False
)

print(f"Saved normalized splicing matrix to: {OUT}")
EOF
