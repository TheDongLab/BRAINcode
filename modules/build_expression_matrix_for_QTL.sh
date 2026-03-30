#!/bin/bash
#SBATCH --job-name=build_expression_matrix
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/build_expression_matrix.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/build_expression_matrix.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0 Python-bundle-PyPI/2024.06-GCCcore-13.3.0

python3 - <<'EOF'
import pandas as pd
from pathlib import Path

BASE = Path("/home/zw529/donglab/data/target_ALS")
OUT  = BASE / "eQTL"
OUT.mkdir(exist_ok=True)

# ─────────────────────────────
# metadata
# ─────────────────────────────
wgs = pd.read_csv(BASE/"targetALS_wgs_metadata.csv")
rna = pd.read_csv(BASE/"targetALS_rnaseq_metadata.csv")

# normalize first (IMPORTANT)
wgs.columns = wgs.columns.str.strip().str.lower()
rna.columns = rna.columns.str.strip().str.lower()

# DEBUG BLOCK
def find_col(df, name):
    return [c for c in df.columns if name in c]

print("WGS subject candidates:", find_col(wgs, "subject"))
print("RNA subject candidates:", find_col(rna, "subject"))

# ─────────────────────────────
# find all normalization files
# ─────────────────────────────
all_tabs = list(BASE.glob("**/RNAseq/Processed/*/normalization.tab"))
print("Found normalization.tab files:", len(all_tabs))

def find_tab(sample_id):
    sid_norm = str(sample_id).replace("-", "_")
    for t in all_tabs:
        parent = t.parent.name
        if sid_norm == parent or sid_norm in parent or str(sample_id) in parent:
            return t
    return None

rna["tab_path"] = rna["externalsampleid"].apply(find_tab)

rna_found = rna[rna["tab_path"].notna()].copy()

print("RNA samples with normalization.tab:", len(rna_found))

# CRITICAL DIAGNOSTIC
print("coverage ratio:", len(rna_found) / max(len(rna), 1))

# ─────────────────────────────
# build TPM matrix
# ─────────────────────────────
expr = None
used_samples = []

for _, row in rna_found.iterrows():
    sample = row["externalsampleid"]
    tab = pd.read_csv(row["tab_path"], sep="\t")

    # strict format check
    if not {"gene_id", "TPM"}.issubset(tab.columns):
        raise ValueError(f"Missing gene_id/TPM in {row['tab_path']}")

    tab = tab[["gene_id", "TPM"]].copy()

    # IMPORTANT: ensure unique sample naming
    colname = str(sample).strip()

    tab.columns = ["gene_id", colname]

    if expr is None:
        expr = tab
    else:
        expr = expr.merge(tab, on="gene_id", how="outer")

    used_samples.append(colname)

print("samples used in matrix:", len(used_samples))
print("expression matrix shape (raw):", expr.shape)

# fill missing values
expr = expr.fillna(0)

# ─────────────────────────────
# genotype alignment
# ─────────────────────────────
geno_raw = BASE / "eQTL/plink/joint_autosomes_matrixEQTL.raw"

geno_header = pd.read_csv(geno_raw, sep=" ", nrows=1)
geno_samples = list(geno_header.columns[1:])

wgs_map = dict(zip(wgs["externalsubjectid"], wgs["externalsampleid"]))

rna_found["geno_id"] = rna_found["externalsubjectid"].map(wgs_map)
rna_found = rna_found.dropna(subset=["geno_id"])

rename_map = dict(zip(rna_found["externalsampleid"], rna_found["geno_id"]))

expr = expr.set_index("gene_id")
expr = expr.rename(columns=rename_map)

# ─────────────────────────────
# HARD DEBUG BEFORE DROP
# ─────────────────────────────
expr_cols = set(expr.columns.astype(str))
geno_set  = set(geno_samples)

print("expr columns:", len(expr_cols))
print("geno columns:", len(geno_set))
print("intersection:", len(expr_cols & geno_set))

# STOP EARLY IF BROKEN
if len(expr_cols & geno_set) == 0:
    raise RuntimeError("NO MATCH BETWEEN EXPRESSION AND GENOTYPE IDS")

# ─────────────────────────────
# final alignment
# ─────────────────────────────
valid_samples = [s for s in geno_samples if s in expr.columns]
expr = expr[valid_samples]

expr = expr.T.groupby(level=0).mean().T

print("Final expression matrix shape (aligned):", expr.shape)

# ─────────────────────────────
# output
# ─────────────────────────────
expr.to_csv(OUT/"expression_matrix.txt", sep="\t")

rna_found[["externalsampleid","externalsubjectid","tissue"]].to_csv(
    OUT/"expression_sample_metadata.csv",
    index=False
)

print("Saved outputs to:", OUT)

EOF
