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
wgs = pd.read_csv(BASE / "targetALS_wgs_metadata.csv")
rna = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")

wgs.columns = wgs.columns.str.strip().str.lower()
rna.columns = rna.columns.str.strip().str.lower()

print("WGS columns:", wgs.columns.tolist())
print("RNA columns:", rna.columns.tolist())

# ─────────────────────────────
# genotype sample IDs from header only
# ─────────────────────────────
geno_raw = BASE / "eQTL/plink/joint_autosomes_matrixEQTL.raw"

import subprocess
first_line = subprocess.check_output(["head", "-n", "1", str(geno_raw)]).decode().strip()
all_geno_cols = first_line.split()
print(f"Total columns in .raw header: {len(all_geno_cols)}")
print(f"First 10 .raw header cols: {all_geno_cols[:10]}")

# PLINK .raw: FID IID PAT MAT SEX PHENOTYPE <SNPs...>
# IID (index 1) is the sample identifier
geno_iid = [norm(x) for x in all_geno_cols[1:6]]
print(f"Cols 1-5 (should be IID,PAT,MAT,SEX,PHENOTYPE): {geno_iid}")

geno_samples = [norm(x) for x in all_geno_cols[6:]]
print(f"Geno samples (SNP columns, first 10): {geno_samples[:10]}")

# Re-read with IID as the sample axis — PLINK .raw uses IID not FID
# Actually collect IIDs row by row from the file using only cols 0+1
geno_meta = pd.read_csv(geno_raw, sep=r"\s+", usecols=["FID", "IID"], engine="python")
geno_iids = [norm(x) for x in geno_meta["IID"].tolist()]
print(f"Geno IIDs (first 10): {geno_iids[:10]}")
print(f"Total geno IIDs: {len(geno_iids)}")

# ─────────────────────────────
# wgs_map: subjectid -> sampleid  AND  sampleid -> sampleid (both directions)
# ─────────────────────────────
wgs_clean = wgs.dropna(subset=["externalsubjectid", "externalsampleid"]).drop_duplicates("externalsubjectid")

wgs_subject_to_sample = {
    norm(r["externalsubjectid"]): norm(r["externalsampleid"])
    for _, r in wgs_clean.iterrows()
}
wgs_sample_to_sample = {
    norm(r["externalsampleid"]): norm(r["externalsampleid"])
    for _, r in wgs_clean.iterrows()
}

print(f"\nWGS subject->sample map size: {len(wgs_subject_to_sample)}")
print(f"WGS subject IDs (first 5): {list(wgs_subject_to_sample.keys())[:5]}")
print(f"WGS sample IDs (first 5):  {list(wgs_subject_to_sample.values())[:5]}")

# ─────────────────────────────
# find normalization files
# ─────────────────────────────
all_tabs = list(BASE.glob("**/RNAseq/Processed/*/normalization.tab"))
print(f"\nFound normalization.tab files: {len(all_tabs)}")

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

print(f"\nSamples used: {len(used_samples)}")
print(f"Expression shape (raw): {expr.shape}")
print(f"Expression sample IDs (first 10): {used_samples[:10]}")

expr = expr.fillna(0).infer_objects(copy=False)

# ─────────────────────────────
# map RNA sample IDs -> genotype IDs
# try subject->sample first, then sample->sample, then direct
# ─────────────────────────────
def resolve_geno_id(rna_subject, rna_sample):
    # path 1: RNA subjectid -> WGS subjectid -> WGS sampleid
    gid = wgs_subject_to_sample.get(norm(rna_subject))
    if gid:
        return gid
    # path 2: RNA sampleid directly in WGS sampleids
    gid = wgs_sample_to_sample.get(norm(rna_sample))
    if gid:
        return gid
    return None

rna_found = rna_found.copy()
rna_found["geno_id"] = rna_found.apply(
    lambda r: resolve_geno_id(r["externalsubjectid"], r["externalsampleid"]), axis=1
)

print(f"\nRNA rows with resolved geno_id: {rna_found['geno_id'].notna().sum()} / {len(rna_found)}")
print("Sample of resolved mapping:")
print(rna_found[["externalsampleid", "externalsubjectid", "geno_id"]].head(10).to_string())

rna_found = rna_found.dropna(subset=["geno_id"])

rename_map = {
    norm(r["externalsampleid"]): norm(r["geno_id"])
    for _, r in rna_found.iterrows()
}

print(f"\nRename map size: {len(rename_map)}")
print(f"Rename map sample (first 5): {list(rename_map.items())[:5]}")

# ─────────────────────────────
# align expression to genotype
# ─────────────────────────────
expr = expr.set_index("gene_id")
expr = expr.rename(columns=rename_map)
expr.columns = [norm(c) for c in expr.columns]

expr_cols = set(expr.columns)
geno_set  = set(geno_iids)

print(f"\nExpr columns:  {len(expr_cols)}  (first 10): {sorted(expr_cols)[:10]}")
print(f"Geno IIDs:     {len(geno_set)}  (first 10): {sorted(geno_set)[:10]}")
print(f"Intersection:  {len(expr_cols & geno_set)}")

if len(expr_cols & geno_set) == 0:
    raise RuntimeError("NO MATCH BETWEEN EXPRESSION AND GENOTYPE IDS")

valid_samples = [s for s in geno_iids if s in expr.columns]
expr = expr[valid_samples]
expr = expr.T.groupby(level=0).mean().T

print(f"Final expression matrix shape: {expr.shape}")

# ─────────────────────────────
# output
# ─────────────────────────────
expr.to_csv(OUT / "expression_matrix.txt", sep="\t")

rna_found[["externalsampleid", "externalsubjectid", "geno_id", "tissue"]].to_csv(
    OUT / "expression_sample_metadata.csv", index=False
)

print(f"Saved outputs to: {OUT}")

EOF
