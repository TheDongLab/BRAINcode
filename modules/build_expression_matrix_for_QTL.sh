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

# ── Load metadata
wgs = pd.read_csv(BASE/"targetALS_wgs_metadata.csv")
rna = pd.read_csv(BASE/"targetALS_rnaseq_metadata.csv")

wgs.columns = wgs.columns.str.strip()
rna.columns = rna.columns.str.strip()
rna.columns = [c.lower() for c in rna.columns]

# ── Shared patients ONLY
shared = set(wgs["Externalsubjectid"].dropna()) & set(rna["externalsubjectid"].dropna())
rna = rna[rna["externalsubjectid"].isin(shared)].copy()

print(f"Patients with BOTH WGS + RNAseq: {len(shared)}")
print(f"RNAseq samples (pre-filter): {len(rna)}")

# ── Find normalization.tab files
all_tabs = list(BASE.glob("**/RNAseq/Processed/*/normalization.tab"))

def find_tab(sample_id):
    sid_norm = sample_id.replace("-", "_")
    for t in all_tabs:
        parent = t.parent.name
        if sid_norm in parent or sample_id in parent:
            return t
    return None

rna["tab_path"] = rna["externalsampleid"].apply(find_tab)

rna_found = rna[rna["tab_path"].notna()].copy()
rna_found = rna_found.drop_duplicates(subset=["externalsampleid"])

print(f"RNAseq samples with normalization.tab: {len(rna_found)}")

# ── Build TPM matrix
expr = None

for _, row in rna_found.iterrows():
    sample = row["externalsampleid"]
    tab = pd.read_csv(row["tab_path"], sep="\t")

    if not {"gene_id","TPM"}.issubset(tab.columns):
        raise ValueError(f"Missing expected columns in {row['tab_path']}")

    tab = tab[["gene_id", "TPM"]].copy()
    tab.columns = ["gene_id", sample]

    if expr is None:
        expr = tab
    else:
        expr = expr.merge(tab, on="gene_id", how="outer")

expr = expr.fillna(0)

print(f"Expression matrix shape: {expr.shape}")

# ── Align sample order with genotype (.raw)
geno_raw = BASE / "eQTL/plink/joint_autosomes_matrixEQTL.raw"

geno_header = pd.read_csv(geno_raw, sep=" ", nrows=1)
geno_samples = list(geno_header.columns[1:])  # drop FID

expr = expr.set_index("gene_id")

common = [s for s in geno_samples if s in expr.columns]
expr = expr[common]

print(f"Final expression matrix shape (aligned): {expr.shape}")

# ── Output
expr_out = OUT / "expression_matrix.txt"
meta_out = OUT / "expression_sample_metadata.csv"

expr.to_csv(expr_out, sep="\t")
rna_found[["externalsampleid","externalsubjectid","tissue"]].to_csv(meta_out, index=False)

print(f"Saved expression matrix: {expr_out}")
print(f"Saved sample metadata: {meta_out}")

EOF
