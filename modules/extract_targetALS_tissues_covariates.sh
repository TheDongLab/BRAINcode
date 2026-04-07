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
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail
module load SAMtools

# ── PATHS ─────────────────────────────────────────────────────────────────────
DATA_DIR="/home/zw529/donglab/data/target_ALS"
OUTDIR="$DATA_DIR/QTL"
mkdir -p "$OUTDIR"

RNAQC_DIR="$OUTDIR/RNAQC_data"
mkdir -p "$RNAQC_DIR"

METADATA="$DATA_DIR/targetALS_rnaseq_metadata.csv"
WGS_META="$DATA_DIR/targetALS_wgs_metadata.csv"
GENCODE_BED="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6"

# Output Files
SKEW_DATA="$RNAQC_DIR/calculated_skew.tsv"
TISSUE_SUMMARY="$OUTDIR/tissue_summary.txt"
PATIENT_TISSUES="$OUTDIR/patient_tissue_breakdown.tsv"
COVARIATE_SUMMARY="$OUTDIR/covariate_summary.txt"
FINAL_COVARIATES="$OUTDIR/covariates.tsv"

# ── STEP 0: CALCULATE RNA SKEWNESS ──────────────────────────────────────────
echo "Running RNA Degradation Rescue (Skewness calculation)..."

TARGET_BED="$RNAQC_DIR/reference_subset.bed"
if [ ! -f "$TARGET_BED" ]; then
    awk -v OFS="\t" '{print $0, $3-$2}' "$GENCODE_BED" | sort -k7,7rn | head -n 500 | cut -f1-6 > "$TARGET_BED"
fi

python3 - <<EOF
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path

def get_skew_for_bam(bam_path, bed_path):
    cmd = f"samtools depth -b {bed_path} {bam_path}"
    try:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
        bins = np.zeros(101)
        genes = []
        with open(bed_path) as f:
            for line in f:
                p = line.split()
                genes.append((p[0], int(p[1]), int(p[2])))
        
        for line in proc.stdout:
            chrom, pos, depth = line.split()
            pos, depth = int(pos), int(depth)
            for g_chr, g_start, g_end in genes:
                if chrom == g_chr and g_start <= pos <= g_end:
                    rel_pos = int(((pos - g_start) / (g_end - g_start)) * 100)
                    if 0 <= rel_pos <= 100: bins[rel_pos] += depth
                    break
        proc.wait()
        if np.sum(bins) == 0: return np.nan
        pos_arr = np.arange(101)
        mean = np.average(pos_arr, weights=bins)
        std = np.sqrt(np.average((pos_arr - mean)**2, weights=bins))
        return np.average(((pos_arr - mean) / std)**3, weights=bins)
    except: return np.nan

if not Path("$SKEW_DATA").exists():
    bams = list(Path("$DATA_DIR").rglob("*.bam"))
    results = []
    for b in bams:
        sid = b.name.split('.')[0].split('_')[0].replace("_", "-")
        results.append({'externalsampleid': sid, 'rna_skew': get_skew_for_bam(str(b), "$TARGET_BED")})
    pd.DataFrame(results).to_csv("$SKEW_DATA", sep="\t", index=False)
EOF

# ── STEP 1: COMPREHENSIVE PROCESSING (TISSUE + COVARIATES) ───────────────────
echo "Generating Tissue Breakdown and Covariate tables..."

python3 - <<EOF
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter

# --- Mappings ---
TISSUE_REMAP = {
    "Motor_Cortex_Lateral": "Motor_Cortex", "Motor_Cortex_Medial": "Motor_Cortex",
    "Medial_Motor_Cortex": "Motor_Cortex", "Lateral_Motor_Cortex": "Motor_Cortex",
    "Lateral_motor_cortex": "Motor_Cortex", "Cortex_Motor_BA4": "Motor_Cortex",
    "BA4_Motor_Cortex": "Motor_Cortex", "Cortex_Motor_Unspecified": "Motor_Cortex",
    "Primary_Motor_Cortex_L": "Motor_Cortex", "Primary_Motor_Cortex_M": "Motor_Cortex",
    "Motor Cortex": "Motor_Cortex", "Motor Cortex Lateral": "Motor_Cortex",
    "Lumbar_spinal_cord": "Lumbar_Spinal_Cord", "Spinal_Cord_Lumbar": "Lumbar_Spinal_Cord",
    "Lumbosacral_Spinal_Cord": "Lumbar_Spinal_Cord", "Lumbar Spinal Cord": "Lumbar_Spinal_Cord",
    "Cervical_spinal_cord": "Cervical_Spinal_Cord", "Spinal_Cord_Cervical": "Cervical_Spinal_Cord",
    "Cervical Spinal Cord": "Cervical_Spinal_Cord", "Frontal Cortex": "Frontal_Cortex"
}

SUBJECT_GROUP_REMAP = {
    "ALS Spectrum MND, Other Neurological Diseases": "ALS Spectrum MND, Other Neurological Disorders",
    "Non Neurological Control": "Non-Neurological Control",
}

C9_REMAP = {
    "ND": "Unknown", "yes": "Yes", "Negative": "No",
    "Not Applicable": "Not Applicable/NaN", "nan": "Not Applicable/NaN"
}

def clean_motor_onset(val):
    val = str(val).strip()
    if val.lower() in ['nan', 'not applicable']: return "Not Applicable/NaN"
    if any(ex in val for ex in ["Bulbar and Limb", "Bulbar/Limb"]): return "Bulbar and Limb"
    if "Axial and Limb" in val: return val
    if "limb" in val.lower() or "extremity" in val.lower(): return "Limb"
    return val

# --- Load ---
df = pd.read_csv("$METADATA")
df.columns = df.columns.str.strip()
wgs_meta = pd.read_csv("$WGS_META")
wgs_meta.columns = wgs_meta.columns.str.strip()

# --- Tissue Tracking ---
wgs_subs = set(wgs_meta["Externalsubjectid"].dropna())
rna_subs = set(df["externalsubjectid"].dropna())
shared = sorted(wgs_subs & rna_subs)

patient_tissues = defaultdict(set)
sample_map = {row['externalsampleid']: row['externalsubjectid'] for _, row in df.iterrows()}
processed_dirs = list(Path("$DATA_DIR").glob("*/RNAseq/Processed/*/"))
for d in processed_dirs:
    raw_tissue = d.parts[-4]
    sid = d.name.replace("_", "-")
    sub_id = sample_map.get(sid)
    if sub_id in shared: patient_tissues[sub_id].add(raw_tissue)

with open("$PATIENT_TISSUES", "w") as f:
    f.write("subject_id\tn_tissues\ttissues\n")
    for sub in sorted(patient_tissues.keys()):
        ts = sorted(patient_tissues[sub])
        f.write(f"{sub}\t{len(ts)}\t{'; '.join(ts)}\n")

all_mapped_tissues = [TISSUE_REMAP.get(t, t.replace(" ", "_")) for ts in patient_tissues.values() for t in ts]
t_counts = Counter(all_mapped_tissues)
with open("$TISSUE_SUMMARY", "w") as f:
    f.write(f"==============================\n target_ALS Tissue Summary\n Shared Patients: {len(shared)}\n==============================\n\n")
    for t, c in sorted(t_counts.items(), key=lambda x: -x[1]):
        f.write(f"{c:<8} {t}\n")

# --- Cleaning ---
df['subject_group'] = df['subject_group'].astype(str).str.strip().replace(SUBJECT_GROUP_REMAP)
df['tissue'] = df['tissue'].map(lambda x: TISSUE_REMAP.get(str(x).strip(), str(x).replace(" ", "_")))
df['c9orf72_repeat_expansion'] = df['c9orf72_repeat_expansion'].astype(str).str.strip().replace(C9_REMAP)
df['site_of_motor_onset'] = df['site_of_motor_onset'].apply(clean_motor_onset)

# Sex Imputation
df['sex'] = df['sex'].astype(str).str.capitalize()
mask = df['sex'].isin(['Unknown', 'Nd', 'Nan', 'None', '']) | df['sex'].isna()
def impute_sex(row):
    gen = str(row.get('sex_genotype', '')).upper()
    return 'Female' if 'XX' in gen else 'Male' if 'XY' in gen else row['sex']
df.loc[mask, 'sex'] = df[mask].apply(impute_sex, axis=1)

# Duration
def calc_dur(row):
    dur = pd.to_numeric(row.get("disease_duration_in_months"), errors='coerce')
    if pd.isna(dur) or dur == 0.0:
        death = pd.to_numeric(row.get("age_at_death"), errors='coerce')
        onset = pd.to_numeric(row.get("age_at_symptom_onset"), errors='coerce')
        if pd.notna(death) and pd.notna(onset):
            return (death - onset) * 12 if death >= onset else np.nan
    return dur
df["disease_duration_in_months"] = df.apply(calc_dur, axis=1)

# --- Rescue Logic ---
if Path("$SKEW_DATA").exists():
    sk_df = pd.read_csv("$SKEW_DATA", sep="\t")
    df = df.merge(sk_df, on='externalsampleid', how='left')

df['rin_num'] = pd.to_numeric(df['rin'], errors='coerce')
df['sk_num']  = pd.to_numeric(df['rna_skew'], errors='coerce')

def qtl_filter(row):
    if pd.notna(row['rin_num']): return True
    if pd.notna(row['sk_num']): return -0.5 <= row['sk_num'] <= 0.5
    return False

df['keep_for_qtl'] = df.apply(qtl_filter, axis=1)

# --- Full Output (Restored ALL missing columns) ---
FULL_COLS = [
    "externalsampleid", "externalsubjectid", "sex", "subject_group", 
    "subject_group_subcategory", "age_at_death", "tissue", "rin", 
    "rna_skew", "keep_for_qtl", "disease_duration_in_months", 
    "post_mortem_interval_in_hours", "site_of_motor_onset", 
    "c9orf72_repeat_expansion", "atxn2_repeat_expansion", "sex_genotype",
    "pct_african", "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas"
]

available = [c for c in FULL_COLS if c in df.columns]
df[available].to_csv("$FINAL_COVARIATES", sep="\t", index=False)

# --- Summary Report ---
with open("$COVARIATE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(f" target_ALS Final Covariate Summary | Pass QTL Filter: {df['keep_for_qtl'].sum()}\n")
    f.write("============================================================\n\n")
    for col in ["sex", "subject_group", "tissue", "keep_for_qtl", "site_of_motor_onset"]:
        f.write(f"── {col.upper()} ──────────────────────────────\n")
        cts = df[col].value_counts(dropna=False)
        for v, c in cts.items():
            f.write(f"    {str(v):<45} {c:>5} ({100*c/len(df):.1f}%)\n")
        f.write("\n")
EOF

echo "Done. All files (Tissue Breakdown, Summary, Full Covariates) generated in $OUTDIR"
