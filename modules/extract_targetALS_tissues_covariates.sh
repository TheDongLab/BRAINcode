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
#   3. Covariate summary (Detailed counts, binned ancestry, and numeric stats)
#   4. Raw covariate extraction per sample
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

# ── Tissue name remapping for directory walking
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
    "Lumbosacaral spinal cord":  "Lumbar_Spinal_Cord",
    "Spinal_Cord_Cervical":      "Cervical_Spinal_Cord",
    "Cervical_spinal_cord":      "Cervical_Spinal_Cord",
    "Spinal_cord_Cervical":      "Cervical_Spinal_Cord",
}

rnaseq_meta = pd.read_csv(RNASEQ_META)
wgs_meta    = pd.read_csv(WGS_META)
rnaseq_meta.columns = rnaseq_meta.columns.str.strip()
wgs_meta.columns    = wgs_meta.columns.str.strip()

wgs_patients = set(wgs_meta["Externalsubjectid"].dropna())
rna_patients = set(rnaseq_meta["externalsubjectid"].dropna())
shared_patients = sorted(wgs_patients & rna_patients)

sample_to_subject = {
    row["externalsampleid"]: row["externalsubjectid"]
    for _, row in rnaseq_meta.iterrows()
    if pd.notna(row["externalsampleid"]) and pd.notna(row["externalsubjectid"])
}

rnaseq_dirs = list(BASE.glob("*/RNAseq/Processed/*/"))
patient_tissues_orig = defaultdict(set)

for d in rnaseq_dirs:
    if not d.is_dir(): continue
    raw_tissue = d.parts[-4]
    sample_id  = d.name.replace("_", "-")
    subject_id = sample_to_subject.get(sample_id)
    if subject_id in shared_patients:
        patient_tissues_orig[subject_id].add(raw_tissue)

with open(OUTDIR / "patient_tissue_breakdown.tsv", "w") as f:
    f.write("subject_id\tn_tissues\ttissues\n")
    for patient in sorted(patient_tissues_orig.keys()):
        tissues = sorted(patient_tissues_orig[patient])
        f.write(f"{patient}\t{len(tissues)}\t{'; '.join(tissues)}\n")

all_tissues = [TISSUE_REMAP.get(t, t.replace(" ", "_")) for tissues in patient_tissues_orig.values() for t in tissues]
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

# ── Mappings for grouping
SUBJECT_GROUP_REMAP = {
    "ALS Spectrum MND, Other Neurological Diseases": "ALS Spectrum MND, Other Neurological Disorders",
    "Non Neurological Control": "Non-Neurological Control",
}

TISSUE_REMAP = {
    "Motor Cortex Lateral": "Motor_Cortex", "Motor Cortex Medial": "Motor_Cortex",
    "Medial Motor Cortex": "Motor_Cortex", "Lateral Motor Cortex": "Motor_Cortex",
    "Lateral motor cortex": "Motor_Cortex", "Lateral_motor_cortex": "Motor_Cortex", 
    "Cortex_Motor_BA4": "Motor_Cortex", "BA4 Motor Cortex": "Motor_Cortex", 
    "Cortex_Motor_Unspecified": "Motor_Cortex", "Primary Motor Cortex L": "Motor_Cortex", 
    "Primary Motor Cortex M": "Motor_Cortex", "Motor Cortex": "Motor_Cortex", 
    "Lumbar spinal cord": "Lumbar_Spinal_Cord", "Spinal Cord Lumbar": "Lumbar_Spinal_Cord", 
    "Lumbosacral Spinal Cord": "Lumbar_Spinal_Cord", "Spinal_Cord_Lumbosacral": "Lumbar_Spinal_Cord", 
    "Lumbar Spinal Cord": "Lumbar_Spinal_Cord", "Lumbosacaral spinal cord": "Lumbar_Spinal_Cord",
    "Cervical Spinal Cord": "Cervical_Spinal_Cord", "Spinal Cord Cervical": "Cervical_Spinal_Cord", 
    "Cervical spinal cord": "Cervical_Spinal_Cord", "Spinal cord Cervical": "Cervical_Spinal_Cord", 
    "Frontal Cortex": "Frontal_Cortex", "Cerebellum": "Cerebellum", 
    "Thoracic Spinal Cord": "Thoracic_Spinal_Cord", "Cortex_Occipital": "Cortex_Occipital", 
    "Choroid": "Choroid", "Liver": "Liver", "Cortex_Temporal": "Cortex_Temporal", 
    "Hippocampus": "Hippocampus", "Cortex_Sensory": "Cortex_Sensory", 
    "Medulla": "Medulla", "Cortex_Unspecified": "Cortex_Unspecified"
}

ANCESTRY_BINS   = [0, 0.01, 0.05, 0.25, 0.50, 0.75, 1.01]
ANCESTRY_LABELS = ["<1%", "1–5%", "5–25%", "25–50%", "50–75%", "75–100%"]

COVARIATES = [
    "externalsampleid", "externalsubjectid", "sex", "subject_group",
    "subject_group_subcategory", "age_at_death", "tissue", "rin",
    "disease_duration_in_months", "post_mortem_interval_in_hours",
    "site_of_motor_onset", "c9orf72_repeat_expansion", "atxn2_repeat_expansion",
    "sex_genotype", "pct_african", "pct_south_asian", "pct_east_asian",
    "pct_european", "pct_americas", "revised_el_escorial_criteria",
    "phenotype_near_time_of_death", "reported_genomic_mutations"
]

df = pd.read_csv(METADATA)
df.columns = df.columns.str.strip()

# ── Clean Data
df['subject_group'] = df['subject_group'].astype(str).str.strip().replace(SUBJECT_GROUP_REMAP)
df['tissue'] = df['tissue'].map(lambda x: TISSUE_REMAP.get(str(x).strip(), str(x).replace(" ", "_")))

# Sex Imputation
df['sex'] = df['sex'].astype(str).str.capitalize()
mask_missing_sex = df['sex'].isin(['Unknown', 'Nd', 'Nan', 'None', '']) | df['sex'].isna()
def impute_sex(row):
    gen = str(row.get('sex_genotype', '')).upper()
    return 'Female' if 'XX' in gen else 'Male' if 'XY' in gen else row['sex']
df.loc[mask_missing_sex, 'sex'] = df[mask_missing_sex].apply(impute_sex, axis=1)

# Duration Calculation
def calculate_duration(row):
    dur = pd.to_numeric(row.get("disease_duration_in_months"), errors='coerce')
    if pd.isna(dur) or dur == 0.0:
        death = pd.to_numeric(row.get("age_at_death"), errors='coerce')
        onset = pd.to_numeric(row.get("age_at_symptom_onset"), errors='coerce')
        if pd.notna(death) and pd.notna(onset):
            return (death - onset) * 12 if death >= onset else np.nan
    return dur
df["disease_duration_in_months"] = df.apply(calculate_duration, axis=1)

# Save Raw TSV
available = [c for c in COVARIATES if c in df.columns]
df[available].to_csv(OUTFILE, sep="\t", index=False)

# ── Generate Detailed Summary
with open(SUMMARY, "w") as f:
    f.write("============================================================\n")
    f.write(" target_ALS Covariate Summary\n")
    f.write(f" Total samples: {len(df)} | Total subjects: {df['externalsubjectid'].nunique()}\n")
    f.write("============================================================\n\n")

    # 1. Categorical
    for col in ["sex", "subject_group", "tissue", "site_of_motor_onset", "c9orf72_repeat_expansion"]:
        if col not in df.columns: continue
        f.write(f"── {col.upper()} ──────────────────────────────\n")
        counts = df[col].value_counts(dropna=False)
        for val, cnt in counts.items():
            f.write(f"    {str(val):<45} {cnt:>5} ({100*cnt/len(df):.1f}%)\n")
        f.write("\n")

    # 2. Numerical (Detailed Stats)
    f.write("── NUMERICAL COVARIATES ────────────────────────────────\n\n")
    for col in ["age_at_death", "rin", "disease_duration_in_months", "post_mortem_interval_in_hours"]:
        if col not in df.columns: continue
        s = pd.to_numeric(df[col], errors='coerce').dropna()
        n_missing = len(df) - len(s)
        if len(s) > 0:
            f.write(f"{col} (n={len(s)}, {n_missing} missing/non-numeric)\n")
            f.write(f"    mean   : {s.mean():.2f} | median : {s.median():.2f} | std : {s.std():.2f}\n")
            f.write(f"    min    : {s.min():.2f} | max    : {s.max():.2f}\n")
            f.write(f"    25th % : {s.quantile(0.25):.2f} | 75th % : {s.quantile(0.75):.2f}\n\n")

    # 3. Ancestry (Binned Distribution)
    f.write("── ANCESTRY COVARIATES (proportion 0–1) ────────────────\n\n")
    for col in ["pct_african", "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas"]:
        if col not in df.columns: continue
        s = pd.to_numeric(df[col], errors='coerce').dropna()
        # Filter for valid proportion range to avoid skewed stats from placeholder values
        s_valid = s[(s >= 0) & (s <= 1)]
        n_invalid = len(s) - len(s_valid)
        
        f.write(f"{col} (n={len(s_valid)} valid, {n_invalid} placeholder/out-of-range)\n")
        f.write(f"    mean: {s_valid.mean():.4f} | median: {s_valid.median():.4f}\n")
        
        binned = pd.cut(s_valid, bins=ANCESTRY_BINS, labels=ANCESTRY_LABELS, right=False)
        for label, cnt in binned.value_counts().sort_index().items():
            f.write(f"      {label:<12} {cnt:>5} ({100*cnt/len(s_valid):.1f}%)\n")
        f.write("\n")

print(f"Summary written: {SUMMARY}")
EOF

echo ""
echo "Done."
echo "  Tissue summary       : $TISSUE_SUMMARY"
echo "  Per-patient breakdown: $PATIENT_TISSUES"
echo "  Covariate summary    : $COVARIATE_SUMMARY"
echo "  Covariates (raw)     : $COVARIATES"
