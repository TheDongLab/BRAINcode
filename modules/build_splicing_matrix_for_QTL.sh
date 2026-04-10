#!/bin/bash
#SBATCH --job-name=build_splicing_matrix
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/build_splicing_matrix.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/build_splicing_matrix.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0 Python-bundle-PyPI/2024.06-GCCcore-13.3.0

python3 - <<'EOF'
import pandas as pd
from pathlib import Path
pd.set_option("future.no_silent_downcasting", True)
BASE = Path("/home/zw529/donglab/data/target_ALS")
OUT  = BASE / "sQTL"
OUT.mkdir(exist_ok=True)
def norm(x):
    return str(x).strip().replace("-", "_")
    
# ─────────────────────────────
# metadata
# ─────────────────────────────
rna = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")
rna.columns = rna.columns.str.strip().str.lower()
print(f"RNA samples in metadata: {len(rna)}")

# ─────────────────────────────
# find PSI files
# ─────────────────────────────
all_psi = list(BASE.glob("**/leafcutter/psi/*.leafcutter.PSI.tsv"))
print(f"Found PSI files: {len(all_psi)}")
def find_psi(sample_id):
    sid_norm = norm(sample_id)
    for p in all_psi:
        # extract sample ID from filename: SAMPLE_ID.leafcutter.PSI.tsv
        fname_sample = norm(p.stem.split(".leafcutter")[0])
        if sid_norm == fname_sample or sid_norm in fname_sample:
            return p
    return None
rna["psi_path"] = rna["externalsampleid"].apply(find_psi)
rna_found = rna[rna["psi_path"].notna()].copy()
print(f"RNA samples with PSI files: {len(rna_found)}")
print(f"Coverage ratio: {len(rna_found) / max(len(rna), 1):.2f}")

# ─────────────────────────────
# build splicing matrix
# ─────────────────────────────
expr = None
used_samples = []
for _, row in rna_found.iterrows():
    sample = norm(row["externalsampleid"])
    psi = pd.read_csv(row["psi_path"], sep="\t")
    if "PSI" not in psi.columns:
        raise ValueError(f"Missing PSI column in {row['psi_path']}")
    # create unique junction ID from genomic coords
    psi["junction_id"] = (
        psi["chrom"] + ":" + psi["strand"] + ":" + 
        psi["start"].astype(str) + "-" + psi["end"].astype(str)
    )
    psi = psi[["junction_id", "PSI"]].copy()
    psi.columns = ["junction_id", sample]
    expr = psi if expr is None else expr.merge(psi, on="junction_id", how="outer")
    used_samples.append(sample)
print(f"Samples in splicing matrix: {len(used_samples)}")
print(f"Splicing matrix shape: {expr.shape}")
expr = expr.fillna(0).infer_objects(copy=False)

# ─────────────────────────────
# output
# ─────────────────────────────
expr.to_csv(OUT / "splicing_matrix.txt", sep="\t", index=False)
rna_found[["externalsampleid", "externalsubjectid", "tissue"]].to_csv(
    OUT / "splicing_sample_metadata.csv", index=False
)
print(f"Saved outputs to: {OUT}")
EOF
