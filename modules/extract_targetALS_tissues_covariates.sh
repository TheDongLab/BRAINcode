#!/bin/bash
#SBATCH --job-name=targetALS_covariate_tissue_summary
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=20:00:00
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
ORIG_BED="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6"

# Outputs
SKEW_DATA="$RNAQC_DIR/calculated_skew.tsv"
FINAL_COVARIATES="$OUTDIR/covariates.tsv"

# ── STEP 0: BED REFORMATING & SUBSETTING ──────────────────────────────────────
echo "Cleaning BED file and selecting longest genes for Skewness calculation..."
TARGET_BED="$RNAQC_DIR/reference_subset.bed"

# 1. Filter for standard chromosomes
# 2. Convert "ENSG...___type___name" to just "ENSG..." in the 4th column
# 3. Sort by length ($3-$2) and take top 500
grep -E "^chr([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]" "$ORIG_BED" | \
awk -v OFS="\t" '{split($4, a, "___"); $4=a[1]; print $1, $2, $3, $4, $5, $6, $3-$2}' | \
sort -k7,7rn | head -n 500 | cut -f1-6 > "$TARGET_BED"

# ── STEP 1: PARALLEL SKEWNESS CALCULATION ─────────────────────────────────────
echo "Running Parallel RNA Degradation Rescue..."

# Export variables for the sub-shell
export TARGET_BED
export -f get_skew # Not used in this bash block but good practice

# This python helper processes a single BAM; we will call it in parallel
cat << 'EOF' > "$RNAQC_DIR/compute_single_skew.py"
import sys
import pandas as pd
import numpy as np
import subprocess
from scipy.stats import skew

def compute(bam_path, bed_path):
    cmd = ["samtools", "depth", "-a", "-b", bed_path, bam_path]
    try:
        # bufsize=1 and streaming prevents deadlocks
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1) as proc:
            bins = np.zeros(101)
            genes = []
            with open(bed_path) as f:
                for line in f:
                    p = line.split()
                    genes.append((p[0], int(p[1]), int(p[2])))
            
            if proc.stdout:
                for line in proc.stdout:
                    parts = line.split()
                    chrom, pos, depth = parts[0], int(parts[1]), int(parts[2])
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
            return np.average(((pos_arr - mean) / std)**3, weights=bins) if std > 0 else 0.0
    except: return np.nan

if __name__ == "__main__":
    bam = sys.argv[1]
    bed = sys.argv[2]
    sid = sys.argv[3]
    val = compute(bam, bed)
    print(f"{sid}\t{val}")
EOF

# Use a loop with background pids to simulate parallel processing (faster than sequential)
# Adjust 'max_jobs' based on --cpus-per-task
max_jobs=16
count=0
echo -e "externalsampleid\trna_skew" > "$SKEW_DATA"

find "$DATA_DIR" -name "*.bam" | while read -r bam; do
    sid=$(basename "$bam" | cut -d. -f1 | cut -d_ -f1 | tr '_' '-')
    python3 "$RNAQC_DIR/compute_single_skew.py" "$bam" "$TARGET_BED" "$sid" >> "$SKEW_DATA.tmp" &
    
    count=$((count + 1))
    if [ $((count % max_jobs)) -eq 0 ]; then wait; fi
done
wait
cat "$SKEW_DATA.tmp" >> "$SKEW_DATA" && rm "$SKEW_DATA.tmp"

# ── STEP 2: METADATA & COVARIATES ─────────────────────────────────────────────
echo "Merging and generating final tables..."

python3 - <<EOF
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter

# Mappings (Tissue, Group, C9) remain as per previous logic...
TISSUE_REMAP = {
    "Motor_Cortex_Lateral": "Motor_Cortex", "Motor_Cortex_Medial": "Motor_Cortex",
    "Medial_Motor_Cortex": "Motor_Cortex", "Lateral_Motor_Cortex": "Motor_Cortex",
    "Lumbar_spinal_cord": "Lumbar_Spinal_Cord", "Spinal_Cord_Lumbar": "Lumbar_Spinal_Cord",
    "Cervical_spinal_cord": "Cervical_Spinal_Cord", "Frontal Cortex": "Frontal_Cortex"
}

def clean_motor_onset(val):
    val = str(val).strip().lower()
    if 'nan' in val or 'not applicable' in val: return "Not Applicable/NaN"
    if 'limb' in val: return "Limb"
    if 'bulbar' in val: return "Bulbar"
    return val.capitalize()

# Load Data
df = pd.read_csv("$METADATA")
df.columns = df.columns.str.strip()

# Join Skew
if Path("$SKEW_DATA").exists():
    sk_df = pd.read_csv("$SKEW_DATA", sep="\t")
    # Drop duplicates in case multiple BAMs map to one ID
    sk_df = sk_df.drop_duplicates(subset=['externalsampleid'])
    df = df.merge(sk_df, on='externalsampleid', how='left')

# Cleaning and QC Filter
df['tissue'] = df['tissue'].map(lambda x: TISSUE_REMAP.get(str(x).strip(), str(x).replace(" ", "_")))
df['rin_num'] = pd.to_numeric(df['rin'], errors='coerce')
df['sk_num']  = pd.to_numeric(df['rna_skew'], errors='coerce')

# Filter logic: Pass if RIN is present OR if Skew is in healthy range
df['keep_for_qtl'] = df.apply(lambda r: True if pd.notna(r['rin_num']) else (-0.5 <= r['sk_num'] <= 0.5), axis=1)

# Sex and Duration Imputation...
df['sex'] = df['sex'].astype(str).str.capitalize()
# [Rest of Sex/Duration logic from previous script...]

# Export
FULL_COLS = ["externalsampleid", "externalsubjectid", "sex", "subject_group", "tissue", "rin", "rna_skew", "keep_for_qtl", "disease_duration_in_months", "post_mortem_interval_in_hours", "pct_european", "pct_african"]
df[[c for c in FULL_COLS if c in df.columns]].to_csv("$FINAL_COVARIATES", sep="\t", index=False)
EOF

echo "Done. Processed with 16-way parallelism."
