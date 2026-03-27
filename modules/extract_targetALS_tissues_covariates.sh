#!/bin/bash
#SBATCH --job-name=targetALS_covariate_tissue_summary
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/targetALS_covariate_tissue_summary.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/targetALS_covariate_tissue_summary.err

# ── target_ALS tissue summary + covariate extraction ─────────────────────────
# Script lives in: ~/donglab/pipelines/scripts/eQTL/
# All output goes to: ~/donglab/data/target_ALS/eQTL/
#
# For all patients with both WGS + RNAseq data:
#   1. Per-patient tissue breakdown (which tissues each patient has RNAseq for)
#   2. Global tissue summary (how many samples per tissue across all patients)
#   3. Covariate summary (counts/frequencies for categorical, mean/median for numerical)
#   4. Raw covariate extraction per sample
#
# Tissue remapping (applied to tissue_summary counts only):
#   Motor_Cortex_Lateral, Motor_Cortex_Medial, Medial_Motor_Cortex,
#   Lateral_Motor_Cortex, Lateral_motor_cortex, Cortex_Motor_BA4,
#   BA4_Motor_Cortex, Cortex_Motor_Unspecified, Primary_Motor_Cortex_L,
#   Primary_Motor_Cortex_M  →  Motor_Cortex
#
#   Lumbar_spinal_cord, Lumbosacral_Spinal_Cord  →  Lumbar_Spinal_Cord
#
#   Spinal_Cord_Cervical  →  Cervical_Spinal_Cord
#
# patient_tissue_breakdown.tsv uses original (raw) tissue names.
#
# Usage:
#   sbatch targetALS_covariate_tissue_summary.sh
#   bash targetALS_covariate_tissue_summary.sh   # run interactively
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

DATA_DIR="/home/zw529/donglab/data/target_ALS"
OUTDIR="$DATA_DIR/eQTL"
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
OUTDIR      = BASE / "eQTL"

# ── Tissue name remapping (raw directory name → canonical name)
# Applied to tissue_summary counts only.
# patient_tissue_breakdown.tsv keeps original names.
TISSUE_REMAP = {
    # Motor Cortex variants
    "Motor_Cortex_Lateral":     "Motor_Cortex",
    "Motor_Cortex_Medial":      "Motor_Cortex",
    "Medial_Motor_Cortex":      "Motor_Cortex",
    "Lateral_Motor_Cortex":     "Motor_Cortex",
    "Lateral_motor_cortex":     "Motor_Cortex",
    "Cortex_Motor_BA4":         "Motor_Cortex",
    "BA4_Motor_Cortex":         "Motor_Cortex",
    "Cortex_Motor_Unspecified": "Motor_Cortex",
    "Primary_Motor_Cortex_L":   "Motor_Cortex",
    "Primary_Motor_Cortex_M":   "Motor_Cortex",
    # Lumbar Spinal Cord variants
    "Lumbar_spinal_cord":       "Lumbar_Spinal_Cord",
    "Lumbosacral_Spinal_Cord":  "Lumbar_Spinal_Cord",
    # Cervical Spinal Cord variants
    "Spinal_Cord_Cervical":     "Cervical_Spinal_Cord",
}

# ── Load metadata
rnaseq_meta = pd.read_csv(RNASEQ_META)
wgs_meta    = pd.read_csv(WGS_META)
rnaseq_meta.columns = rnaseq_meta.columns.str.strip()
wgs_meta.columns    = wgs_meta.columns.str.strip()

# ── Identify patients with both WGS + RNAseq
wgs_patients    = set(wgs_meta["Externalsubjectid"].dropna())
rna_patients    = set(rnaseq_meta["externalsubjectid"].dropna())
shared_patients = sorted(wgs_patients & rna_patients)
print(f"Patients with both WGS + RNAseq: {len(shared_patients)}")

# ── Build sample ID -> subject ID lookup (normalize hyphens)
sample_to_subject = {
    row["externalsampleid"]: row["externalsubjectid"]
    for _, row in rnaseq_meta.iterrows()
    if pd.notna(row["externalsampleid"]) and pd.notna(row["externalsubjectid"])
}

# ── Walk RNAseq/Processed/ directories and map to subject + tissue
# Structure: target_ALS / <TISSUE> / RNAseq / Processed / <SAMPLE> /
#                          parts[-4]                        parts[-1]
#
# Two separate accumulators:
#   patient_tissues_raw  – list of canonical names per patient
#                          (allows duplicates so all subtypes are counted)
#   patient_tissues_orig – set of original names per patient
#                          (deduplicated, used for breakdown TSV)
rnaseq_dirs          = list(BASE.glob("*/RNAseq/Processed/*/"))
patient_tissues_raw  = defaultdict(list)
patient_tissues_orig = defaultdict(set)
unmatched            = []

for d in rnaseq_dirs:
    if not d.is_dir():
        continue
    raw_tissue = d.parts[-4]
    canonical  = TISSUE_REMAP.get(raw_tissue, raw_tissue)
    dir_name   = d.name
    sample_id  = dir_name.replace("_", "-")

    subject_id = sample_to_subject.get(sample_id)
    if subject_id is None:
        unmatched.append(dir_name)
        continue

    if subject_id in shared_patients:
        patient_tissues_raw[subject_id].append(canonical)   # list → keeps per-subtype dupes
        patient_tissues_orig[subject_id].add(raw_tissue)    # set  → original names, deduped

print(f"Unmatched directories (not in metadata): {len(unmatched)}")
if unmatched:
    for u in unmatched[:10]:
        print(f"  {u}")
    if len(unmatched) > 10:
        print(f"  ... and {len(unmatched) - 10} more")

# ── Write per-patient tissue breakdown TSV (original names)
with open(OUTDIR / "patient_tissue_breakdown.tsv", "w") as f:
    f.write("subject_id\tn_tissues\ttissues\n")
    for patient in sorted(patient_tissues_orig.keys()):
        tissues = sorted(patient_tissues_orig[patient])
        f.write(f"{patient}\t{len(tissues)}\t{'; '.join(tissues)}\n")

print(f"Per-patient tissue breakdown written: {len(patient_tissues_orig)} patients")

# ── Write global tissue summary (canonical names, all subtypes counted)
all_tissues   = [t for tissues in patient_tissues_raw.values() for t in tissues]
tissue_counts = Counter(all_tissues)

with open(OUTDIR / "tissue_summary.txt", "w") as f:
    f.write("==============================\n")
    f.write(" target_ALS Tissue Summary\n")
    f.write(f" Patients with WGS + RNAseq: {len(shared_patients)}\n")
    f.write(f" Patients with tissue data:  {len(patient_tissues_raw)}\n")
    f.write("==============================\n\n")
    f.write(f"{'COUNT':<8} {'TISSUE'}\n")
    f.write(f"{'-----':<8} {'------------------------------'}\n")
    for tissue, count in sorted(tissue_counts.items(), key=lambda x: -x[1]):
        f.write(f"{count:<8} {tissue}\n")
    f.write(f"\nTotal tissue-sample entries: {sum(tissue_counts.values())}\n")
    f.write(f"Unique tissues:              {len(tissue_counts)}\n")

print(f"Global tissue summary written: {len(tissue_counts)} unique tissues")
EOF

echo ""
echo "Tissue summary       : $TISSUE_SUMMARY"
echo "Per-patient breakdown: $PATIENT_TISSUES"

# ── Step 4: Covariate summary + raw extraction ────────────────────────────────
echo ""
echo "Extracting covariates and generating summary..."

python3 - <<'EOF'
import pandas as pd

DATA_DIR = "/home/zw529/donglab/data/target_ALS"
METADATA  = f"{DATA_DIR}/targetALS_rnaseq_metadata.csv"
OUTFILE   = f"{DATA_DIR}/eQTL/covariates.tsv"
SUMMARY   = f"{DATA_DIR}/eQTL/covariate_summary.txt"

COVARIATES = [
    "externalsampleid",
    "externalsubjectid",
    "sex",
    "ethnicity",
    "subject_group",
    "subject_group_subcategory",
    "age_at_death",
    "age_at_symptom_onset",
    "tissue",
    "rin",
    "ph",
    "disease_duration_in_months",
    "post_mortem_interval_in_hours",
    "site_of_motor_onset",
    "c9orf72_repeat_expansion",
    "atxn2_repeat_expansion",
    "c9_repeat_size",
    "atxn2_repeat_size",
    "sex_genotype",
    "pct_african",
    "pct_south_asian",
    "pct_east_asian",
    "pct_european",
    "pct_americas",
    "revised_el_escorial_criteria",
    "phenotype_near_time_of_death",
    "reported_genomic_mutations",
    "platform",
    "site_specimen_collected",
]

# Shown as value counts + percentages
CATEGORICAL = [
    "sex", "ethnicity", "subject_group", "subject_group_subcategory",
    "tissue", "site_of_motor_onset", "c9orf72_repeat_expansion",
    "atxn2_repeat_expansion", "sex_genotype", "revised_el_escorial_criteria",
    "phenotype_near_time_of_death", "reported_genomic_mutations",
]

# Shown as mean/median/std/min/max/IQR
NUMERICAL = [
    "age_at_death", "age_at_symptom_onset", "rin",
    "post_mortem_interval_in_hours",
]

# Shown as mean/median/std/min/max/IQR + binned distribution (×100 for %)
ANCESTRY = [
    "pct_african", "pct_south_asian", "pct_east_asian",
    "pct_european", "pct_americas",
]

ANCESTRY_BINS   = [0, 0.01, 0.05, 0.25, 0.50, 0.75, 1.01]
ANCESTRY_LABELS = ["<1%", "1–5%", "5–25%", "25–50%", "50–75%", "75–100%"]

# ── Raw covariate extraction
df = pd.read_csv(METADATA)
df.columns = df.columns.str.strip()

available = [c for c in COVARIATES if c in df.columns]
missing   = [c for c in COVARIATES if c not in df.columns]

df[available].to_csv(OUTFILE, sep="\t", index=False)
print(f"Covariates extracted: {len(df)} rows x {len(available)} columns -> {OUTFILE}")
if missing:
    print(f"WARNING: columns not found: {missing}")

# ── Covariate summary report
with open(SUMMARY, "w") as f:

    f.write("=" * 60 + "\n")
    f.write(" target_ALS Covariate Summary\n")
    f.write(f" Total samples: {len(df)}\n")
    f.write(f" Total subjects: {df['externalsubjectid'].nunique()}\n")
    f.write("=" * 60 + "\n\n")

    # Categorical covariates
    f.write("── CATEGORICAL COVARIATES ──────────────────────────────\n\n")
    for col in CATEGORICAL:
        if col not in df.columns:
            continue
        counts    = df[col].value_counts(dropna=False)
        n_missing = df[col].isna().sum()
        f.write(f"{col}  (n={len(df) - n_missing} non-null, {n_missing} missing)\n")
        for val, cnt in counts.items():
            pct = 100 * cnt / len(df)
            f.write(f"    {str(val):<45} {cnt:>5}  ({pct:.1f}%)\n")
        f.write("\n")

    # Numerical covariates
    f.write("── NUMERICAL COVARIATES ────────────────────────────────\n\n")
    for col in NUMERICAL:
        if col not in df.columns:
            continue
        s         = pd.to_numeric(df[col], errors='coerce').dropna()
        n_missing = len(df) - len(s)
        if len(s) == 0:
            f.write(f"{col}  (all missing)\n\n")
            continue
        f.write(f"{col}  (n={len(s)}, {n_missing} missing/non-numeric)\n")
        f.write(f"    mean   : {s.mean():.2f}\n")
        f.write(f"    median : {s.median():.2f}\n")
        f.write(f"    std    : {s.std():.2f}\n")
        f.write(f"    min    : {s.min():.2f}\n")
        f.write(f"    max    : {s.max():.2f}\n")
        f.write(f"    25th % : {s.quantile(0.25):.2f}\n")
        f.write(f"    75th % : {s.quantile(0.75):.2f}\n")
        f.write("\n")

    # Ancestry covariates — stats + binned distribution
    f.write("── ANCESTRY COVARIATES (proportion 0–1) ────────────────\n\n")
    for col in ANCESTRY:
        if col not in df.columns:
            continue
        s         = pd.to_numeric(df[col], errors='coerce').dropna()
        n_missing = len(df) - len(s)
        if len(s) == 0:
            f.write(f"{col}  (all missing)\n\n")
            continue
        f.write(f"{col}  (n={len(s)}, {n_missing} missing/non-numeric)\n")
        f.write(f"    Distribution:\n")
        binned = pd.cut(s, bins=ANCESTRY_BINS, labels=ANCESTRY_LABELS, right=False)
        for label, cnt in binned.value_counts().sort_index().items():
            pct = 100 * cnt / len(s)
            f.write(f"      {label:<12} {cnt:>5}  ({pct:.1f}%)\n")
        f.write("\n")

print(f"Covariate summary written: {SUMMARY}")
EOF

echo ""
echo "Done."
echo "  Tissue summary       : $TISSUE_SUMMARY"
echo "  Per-patient breakdown: $PATIENT_TISSUES"
echo "  Covariate summary    : $COVARIATE_SUMMARY"
echo "  Covariates (raw)     : $COVARIATES"
