#!/bin/bash
#SBATCH --job-name=build_expression_matrix
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/build_expression_matrix.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/build_expression_matrix.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

# Set to target_ALS or answer_ALS before running
EXPERIMENT="answer_ALS"
export EXPERIMENT

unset PYTHONPATH

module purge
module load Python/3.12.3-GCCcore-13.3.0
module load Python-bundle-PyPI/2024.06-GCCcore-13.3.0

python3 - <<'EOF'
import os
import pandas as pd
import numpy as np
from scipy.stats import norm as scipy_norm
from pathlib import Path

pd.set_option("future.no_silent_downcasting", True)

EXPERIMENT = os.environ["EXPERIMENT"]
BASE = Path("/home/zw529/donglab/data") / EXPERIMENT
OUT = BASE / "QTL"
OUT.mkdir(exist_ok=True)

def norm_id(x):
    return str(x).strip().replace("-", "_")

if EXPERIMENT == "target_ALS":
    # Keep the existing TargetALS metadata matching and output columns.
    rna = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")
    rna.columns = rna.columns.str.strip().str.lower()
    print(f"RNA samples in metadata: {len(rna)}")

    all_tabs = list(BASE.glob("**/RNAseq/Processed/*/normalization.tab"))
    print(f"Found normalization.tab files: {len(all_tabs)}")

    def find_tab(sample_id):
        sid_norm = norm_id(sample_id)
        for tab_path in all_tabs:
            parent = norm_id(tab_path.parent.name)
            if sid_norm == parent or sid_norm in parent:
                return tab_path
        return None

    rna["tab_path"] = rna["externalsampleid"].apply(find_tab)
    rna_found = rna[rna["tab_path"].notna()].copy()
    print(f"RNA samples with normalization.tab: {len(rna_found)}")

elif EXPERIMENT == "answer_ALS":
    # AnswerALS processed folders are named by donor and contain normalization.tab.
    rna_processed = BASE / "iPSC_Motor_Neuron/RNAseq/Processed"
    all_tabs = sorted(rna_processed.glob("*/normalization.tab"))
    print(f"Found AnswerALS normalization.tab files: {len(all_tabs)}")

    rna_found = pd.DataFrame([
        {
            "matrix_sample_id": tab.parent.name,
            "donor_id": tab.parent.name,
            "tissue": "iPSC_Motor_Neuron",
            "tab_path": tab,
        }
        for tab in all_tabs
    ])
    print(f"AnswerALS donors with normalization.tab: {len(rna_found)}")

else:
    raise ValueError("EXPERIMENT must be target_ALS or answer_ALS")

# Build the TPM expression matrix.
expr = None
for _, row in rna_found.iterrows():
    if EXPERIMENT == "target_ALS":
        sample = norm_id(row["externalsampleid"])
    else:
        # Preserve CASE-... IDs to match the AnswerALS joint VCF.
        sample = str(row["matrix_sample_id"]).strip()

    tab = pd.read_csv(row["tab_path"], sep="\t")

    if not {"gene_id", "TPM"}.issubset(tab.columns):
        continue

    tab = tab[["gene_id", "TPM"]].copy()
    tab.columns = ["gene_id", sample]

    if expr is None:
        expr = tab
    else:
        expr = expr.merge(tab, on="gene_id", how="outer")

if expr is None:
    raise RuntimeError(f"No usable TPM tables found for {EXPERIMENT}")

expr = expr.fillna(0).infer_objects(copy=False)
print(f"Raw matrix shape: {expr.shape}")

# Rank-based inverse normal transformation across samples for each gene.
def inverse_normal_transform(series):
    n = len(series)
    ranks = series.rank(method="average")
    prob = (ranks - 0.5) / n
    return scipy_norm.ppf(prob)

print("Applying rank-based inverse normal transformation...")
gene_ids = expr["gene_id"]
numeric_data = expr.drop("gene_id", axis=1)
sample_names = numeric_data.columns

transformed = numeric_data.apply(
    lambda x: inverse_normal_transform(x) if x.std() > 1e-9 else np.zeros(len(x)),
    axis=1,
    result_type="expand",
)
transformed.columns = sample_names
expr_final = pd.concat([gene_ids, transformed], axis=1)

expr_final.to_csv(OUT / "expression_matrix.txt", sep="\t", index=False)

if EXPERIMENT == "target_ALS":
    rna_found[["externalsampleid", "externalsubjectid", "tissue"]].to_csv(
        OUT / "expression_sample_metadata.csv", index=False
    )
else:
    rna_found[["matrix_sample_id", "donor_id", "tissue"]].to_csv(
        OUT / "expression_sample_metadata.csv", index=False
    )

print(f"Saved normalized outputs to: {OUT}")
EOF
