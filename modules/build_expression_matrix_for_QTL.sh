#!/bin/bash
#SBATCH --job-name=build_expression_matrix
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/build_expression_matrix.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/build_expression_matrix.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0 Python-bundle-PyPI/2024.06-GCCcore-13.3.0

python3 - <<'EOF'
import pandas as pd
from pathlib import Path

pd.set_option("future.no_silent_downcasting", True)

BASE = Path("/home/zw529/donglab/data/target_ALS")
OUT  = BASE / "eQTL"
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
# find normalization files
# ─────────────────────────────
all_tabs = list(BASE.glob("**/RNAseq/Processed/*/normalization.tab"))
print(f"Found normalization.tab files: {len(all_tabs)}")

def find_tab(sample_id):
    sid_norm = norm(sample_id)
    for t in all_tabs:
        parent = norm(t.parent.name)
        if sid_norm == parent or sid_norm in parent:
            return t
    return None

rna["tab_path"] = rna["externalsampleid"].apply(find_tab)
rna_found = rna[rna["tab_path"].notna()].copy()

print(f"RNA samples with normalization.tab: {len(rna_found)}")
print(f"Coverage ratio: {len(rna_found) / max(len(rna), 1):.2f}")

# ─────────────────────────────
# build expression matrix
# ─────────────────────────────
expr = None
used_samples = []

for _, row in rna_found.iterrows():
    sample = norm(row["externalsampleid"])
    tab = pd.read_csv(row["tab_path"], sep="\t")

    if not {"gene_id", "TPM"}.issubset(tab.columns):
        raise ValueError(f"Missing gene_id/TPM in {row['tab_path']}")

    tab = tab[["gene_id", "TPM"]].copy()
    tab.columns = ["gene_id", sample]

    expr = tab if expr is None else expr.merge(tab, on="gene_id", how="outer")
    used_samples.append(sample)

print(f"Samples in expression matrix: {len(used_samples)}")
print(f"Expression matrix shape: {expr.shape}")

expr = expr.fillna(0).infer_objects(copy=False)

# ─────────────────────────────
# output
# ─────────────────────────────
expr.to_csv(OUT / "expression_matrix.txt", sep="\t", index=False)

rna_found[["externalsampleid", "externalsubjectid", "tissue"]].to_csv(
    OUT / "expression_sample_metadata.csv", index=False
)

print(f"Saved outputs to: {OUT}")

EOF
