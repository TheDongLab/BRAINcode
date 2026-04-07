#!/bin/bash
#SBATCH --job-name=targetALS_covariate_tissue_summary
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/targetALS_covariate_tissue_summary.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/targetALS_covariate_tissue_summary.err

# ── target_ALS tissue summary + covariate extraction ─────────────────────────
# Script lives in: ~/donglab/pipelines/scripts/QTL/
# All output goes to: ~/donglab/data/target_ALS/QTL/
#
# For all patients with both WGS + RNAseq data:
#   1. Per-patient tissue breakdown (which tissues each patient has RNAseq for)
#   2. Global tissue summary (how many samples per tissue across all patients)
#   3. Covariate summary (counts/frequencies for categorical, mean/median for numerical)
#   4. Raw covariate extraction per sample
#
# Tissue remapping (canonical names used in tissue_summary only):
#   Motor_Cortex_Lateral, Motor_Cortex_Medial, Medial_Motor_Cortex,
#   Lateral_Motor_Cortex, Lateral_motor_cortex, Cortex_Motor_BA4,
#   BA4_Motor_Cortex, Cortex_Motor_Unspecified, Primary_Motor_Cortex_L,
#   Primary_Motor_Cortex_M  →  Motor_Cortex
#
#   Lumbar_spinal_cord, Lumbosacral_Spinal_Cord,
#   Spinal_Cord_Lumbar, Spinal_Cord_Lumbosacral  →  Lumbar_Spinal_Cord
#
#   Spinal_Cord_Cervical, Cervical_spinal_cord,
#   Spinal_cord_Cervical                         →  Cervical_Spinal_Cord
#
# patient_tissue_breakdown.tsv uses original (raw) tissue names.
# tissue_summary.txt counts are derived from the same set, guaranteeing
# both outputs sum to the same total.
#
# Usage:
#   sbatch targetALS_covariate_tissue_summary.sh
#   bash targetALS_covariate_tissue_summary.sh   # run interactively
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

DATA_DIR="/home/zw529/donglab/data/target_ALS"
OUTDIR="$DATA_DIR/QTL"
mkdir -p "$OUTDIR"

TISSUE_SUMMARY="$OUTDIR/tissue_summary.txt"
PATIENT_TISSUES="$OUTDIR/patient_tissue_breakdown.tsv"
COVARIATE_SUMMARY="$OUTDIR/covariate_summary.txt"
COVARIATES="$OUTDIR/covariates.tsv"

# ── Steps 1-3: Per-patient tissue breakdown + global summary ──────────────────
echo "Building per-patient tissue breakdown and global summary..."

python3 - <<'EOF'
import pandas as pd
from pathlib import Path
from collections import defaultdict, Counter

BASE        = Path("/home/zw529/donglab/data/target_ALS")
RNASEQ_META = BASE / "targetALS_rnaseq_metadata.csv"
WGS_META    = BASE / "targetALS_wgs_metadata.csv"
OUTDIR      = BASE / "QTL"

# ── Tissue name remapping (raw directory name → canonical name)
TISSUE_REMAP = {
    "Motor_Cortex_Lateral":     "Motor_Cortex",
    "Motor_Cortex_Medial":      "Motor_Cortex",
    "Medial_Motor_Cortex":      "Motor_Cortex",
    "Lateral_Motor_Cortex":     "Motor_Cortex",
    "Lateral_motor_cortex":     "Motor_Cortex",
    "Cortex_Motor_BA4":          "Motor_Cortex",
    "BA4_Motor_Cortex":          "Motor_Cortex",
    "Cortex_Motor_Unspecified": "Motor_Cortex",
    "Primary_Motor_Cortex_L":   "Motor_Cortex",
    "Primary_Motor_Cortex_M":   "Motor_Cortex",
    "Lumbar_spinal_cord":        "Lumbar_Spinal_Cord",
    "Lumbosacral_Spinal_Cord":  "Lumbar_Spinal_Cord",
    "Spinal_Cord_Lumbar":        "Lumbar_Spinal_Cord",
    "Spinal_Cord_Lumbosacral":  "Lumbar_Spinal_Cord",
    "Spinal_Cord_Cervical":      "Cervical_Spinal_Cord",
    "Cervical_spinal_cord":      "Cervical_Spinal_Cord",
    "Spinal_cord_Cervical":      "Cervical_Spinal_Cord",
}

# ── Load metadata
rnaseq_meta = pd.read_csv(RNASEQ_META)
wgs_meta    = pd.read_csv(WGS_META)
rnaseq_meta.columns = rnaseq_meta.columns.str.strip()
wgs_meta.columns    = wgs_meta.columns.str.strip()

wgs_patients    = set(wgs_meta["Externalsubjectid"].dropna())
rna_patients    = set(rnaseq_meta["externalsubjectid"].dropna())
shared_patients = sorted(wgs_patients & rna_patients)

sample_to_subject = {
    row["externalsampleid"]: row["externalsubjectid"]
    for _, row in rnaseq_meta.iterrows()
    if pd.notna(row["externalsampleid"]) and pd.notna(row["externalsubjectid"])
}

rnaseq_dirs          = list(BASE.glob("*/RNAseq/Processed/*/"))
patient_tissues_orig = defaultdict(set)

for d in rnaseq_dirs:
    if not d.is_dir():
        continue
    raw_tissue = d.parts[-4]
    sample_id  = d.name.replace("_", "-")
    subject_id = sample_to_subject.get(sample_id)

    if subject_id in shared_patients:
        patient_tissues_orig[subject_id].add(raw_tissue)

# ── Write per-patient tissue breakdown
with open(OUTDIR / "patient_tissue_breakdown.tsv", "w") as f:
    f.write("subject_id\tn_tissues\ttissues\n")
    for patient in sorted(patient_tissues_orig.keys()):
        tissues = sorted(patient_tissues_orig[patient])
        f.write(f"{patient}\t{len(tissues)}\t{'; '.join(tissues)}\n")

# ── Derive counts
all_tissues   = [TISSUE_REMAP.get(t, t) for tissues in patient_tissues_orig.values() for t in tissues]
tissue_counts = Counter(all_tissues)

with open(OUTDIR / "tissue_summary.txt", "w") as f:
    f.write("==============================\n")
    f.write(" target_ALS Tissue Summary\n")
    f.write(f" Patients with WGS + RNAseq: {len(shared_patients)}\n")
    f.write("==============================\n\n")
    f.write(f"{'COUNT':<8} {'TISSUE'}\n")
    f.write(f"{'-----':<8} {'------------------------------'}\n")
    for tissue, count in sorted(tissue_counts.items(), key=lambda x: -x[1]):
        f.write(f"{count:<8} {tissue}\n")
EOF

# ── Step 4: Covariate summary + raw extraction ────────────────────────────────
echo "Extracting covariates and generating summary..."

python3 - <<'EOF'
import pandas as pd
import numpy as np

DATA_DIR = "/home/zw529/donglab/data/target_ALS"
METADATA  = f"{DATA_DIR}/targetALS_rnaseq_metadata.csv"
OUTFILE   = f"{DATA_DIR}/QTL/covariates.tsv"
SUMMARY   = f"{DATA_DIR}/QTL/covariate_summary.txt"

COVARIATES = [
    "externalsampleid", "externalsubjectid", "sex", "subject_group",
    "subject_group_subcategory", "age_at_death", "tissue", "rin",
    "disease_duration_in_months", "post_mortem_interval_in_hours",
    "site_of_motor_onset", "c9orf72_repeat_expansion", "atxn2_repeat_expansion",
    "c9_repeat_size", "atxn2_repeat_size", "sex_genotype", "pct_african",
    "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas",
    "revised_el_escorial_criteria", "phenotype_near_time_of_death",
    "reported_genomic_mutations", "platform", "site_specimen_collected",
]

CATEGORICAL = [
    "sex", "subject_group", "subject_group_subcategory", "tissue", 
    "site_of_motor_onset", "c9orf72_repeat_expansion", "atxn2_repeat_expansion", 
    "sex_genotype", "revised_el_escorial_criteria", "phenotype_near_time_of_death", 
    "reported_genomic_mutations",
]

NUMERICAL = [
    "age_at_death", "rin", "post_mortem_interval_in_hours", "disease_duration_in_months",
]

ANCESTRY = ["pct_african", "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas"]
ANCESTRY_BINS   = [0, 0.01, 0.05, 0.25, 0.50, 0.75, 1.01]
ANCESTRY_LABELS = ["<1%", "1–5%", "5–25%", "25–50%", "50–75%", "75–100%"]

df = pd.read_csv(METADATA)
df.columns = df.columns.str.strip()

# ── 1. Sex Normalization and Imputation
df['sex'] = df['sex'].astype(str).str.capitalize()
missing_tags = ['Unknown', 'Nd', 'Nan', 'None', '']
mask_missing_sex = df['sex'].isin(missing_tags) | df['sex'].isna()

def impute_sex(row):
    # Ensure genotype is treated as string to avoid errors on XX/XY check
    genotype = str(row.get('sex_genotype', '')).upper()
    if 'XX' in genotype:
        return 'Female'
    elif 'XY' in genotype:
        return 'Male'
    return row['sex']

df.loc[mask_missing_sex, 'sex'] = df[mask_missing_sex].apply(impute_sex, axis=1)

# ── 2. Disease Duration Logic with Type Safety
def calculate_duration(row):
    # Convert duration to numeric first to handle 0.0 correctly
    dur = pd.to_numeric(row.get("disease_duration_in_months"), errors='coerce')
    
    # Check if duration is 0, NaN, or missing
    if pd.isna(dur) or dur == 0.0:
        # Convert age inputs to numeric to avoid TypeErrors during comparison
        death = pd.to_numeric(row.get("age_at_death"), errors='coerce')
        onset = pd.to_numeric(row.get("age_at_symptom_onset"), errors='coerce')
        
        if pd.notna(death) and pd.notna(onset):
            if onset > death:
                return np.nan
            return (death - onset) * 12
    return dur

df["disease_duration_in_months"] = df.apply(calculate_duration, axis=1)

# Save TSV
available = [c for c in COVARIATES if c in df.columns]
df[available].to_csv(OUTFILE, sep="\t", index=False)

# ── Covariate summary report
with open(SUMMARY, "w") as f:
    f.write("============================================================\n")
    f.write(" target_ALS Covariate Summary\n")
    f.write(f" Total samples: {len(df)}\n")
    f.write(f" Total subjects: {df['externalsubjectid'].nunique()}\n")
    f.write("============================================================\n\n")

    f.write("── CATEGORICAL COVARIATES ──────────────────────────────\n\n")
    for col in CATEGORICAL:
        if col not in df.columns: continue
        counts = df[col].value_counts(dropna=False)
        n_missing = df[col].isna().sum()
        f.write(f"{col} (n={len(df) - n_missing} non-null, {n_missing} missing)\n")
        for val, cnt in counts.items():
            f.write(f"    {str(val):<45} {cnt:>5} ({100*cnt/len(df):.1f}%)\n")
        f.write("\n")

    f.write("── NUMERICAL COVARIATES ────────────────────────────────\n\n")
    for col in NUMERICAL:
        if col not in df.columns: continue
        s = pd.to_numeric(df[col], errors='coerce').dropna()
        n_missing = len(df) - len(s)
        f.write(f"{col} (n={len(s)}, {n_missing} missing/non-numeric)\n")
        f.write(f"    mean   : {s.mean():.2f}\n")
        f.write(f"    median : {s.median():.2f}\n")
        f.write(f"    std    : {s.std():.2f}\n")
        f.write(f"    min    : {s.min():.2f}\n")
        f.write(f"    max    : {s.max():.2f}\n\n")

    f.write("── ANCESTRY COVARIATES (proportion 0–1) ────────────────\n\n")
    for col in ANCESTRY:
        if col not in df.columns: continue
        s = pd.to_numeric(df[col], errors='coerce').dropna()
        n_missing = len(df) - len(s)
        f.write(f"{col} (n={len(s)}, {n_missing} missing/non-numeric)\n")
        binned = pd.cut(s, bins=ANCESTRY_BINS, labels=ANCESTRY_LABELS, right=False)
        for label, cnt in binned.value_counts().sort_index().items():
            f.write(f"      {label:<12} {cnt:>5} ({100*cnt/len(s):.1f}%)\n")
        f.write("\n")

print(f"Summary written: {SUMMARY}")
EOF

echo ""
echo "Done."
echo "  Tissue summary       : $TISSUE_SUMMARY"
echo "  Per-patient breakdown: $PATIENT_TISSUES"
echo "  Covariate summary    : $COVARIATE_SUMMARY"
echo "  Covariates (raw)     : $COVARIATES"
