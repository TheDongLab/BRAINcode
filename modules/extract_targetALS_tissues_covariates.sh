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

# Output Files
TARGET_BED="/home/zw529/donglab/data/target_ALS/QTL/RNAQC_data/reference_subset.bed"
SKEW_DATA="$RNAQC_DIR/calculated_skew.tsv"
FINAL_COVARIATES="$OUTDIR/covariates.tsv"
PATIENT_TISSUES="$OUTDIR/patient_tissue_breakdown.tsv"
TISSUE_SUMMARY="$OUTDIR/tissue_summary.txt"
COVARIATE_SUMMARY="$OUTDIR/covariate_summary.txt"

# ── STEP 0: BED REFORMATTING ──────────────────────────────────────────────────
echo 'Cleaning BED file and selecting longest genes...'

# Temporarily disable pipefail if you use 'head' to prevent SIGPIPE crashes
set +o pipefail

grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]' "$ORIG_BED" | \
awk -v OFS='\t' '{split($4, a, "___"); $4=a[1]; print $1, $2, $3, $4, $5, $6, $3-$2}' | \
sort -k7,7rn | \
head -n 500 | \
cut -f1-6 > "$TARGET_BED" || true   # The '|| true' ensures grep failing doesn't kill the script

set -o pipefail # Turn it back on for the rest of the script

echo "Proceeding to next step..."

# ── STEP 1: PARALLEL SKEWNESS CALCULATION ─────────────────────────────────────
echo "Running Parallel RNA Degradation Rescue..."

# Worker script - strictly defined
cat << 'EOF' > "$RNAQC_DIR/worker.py"
import sys, subprocess, numpy as np
def compute(bam, bed):
    cmd = ["samtools", "depth", "-a", "-b", bed, bam]
    try:
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1) as proc:
            bins = np.zeros(101)
            genes = []
            with open(bed) as f:
                for line in f:
                    p = line.split(); genes.append((p[0], int(p[1]), int(p[2])))
            if proc.stdout:
                for line in proc.stdout:
                    parts = line.split()
                    if len(parts) < 3: continue
                    chrom, pos, depth = parts[0], int(parts[1]), int(parts[2])
                    for g_chr, g_start, g_end in genes:
                        if chrom == g_chr and g_start <= pos <= g_end:
                            rel = int(((pos - g_start) / (g_end - g_start)) * 100)
                            if 0 <= rel <= 100: bins[rel] += depth
                            break
            proc.wait()
            if np.sum(bins) == 0: return "NaN"
            pos_arr = np.arange(101)
            mean = np.average(pos_arr, weights=bins)
            std = np.sqrt(np.average((pos_arr - mean)**2, weights=bins))
            return str(np.average(((pos_arr - mean) / std)**3, weights=bins)) if std > 0 else "0.0"
    except: return "NaN"

if __name__ == "__main__":
    bam_path = sys.argv[1]
    bed_path = sys.argv[2]
    sid = sys.argv[3]
    print(f"{sid}\t{compute(bam_path, bed_path)}")
EOF

# Create sample list for xargs
SAMPLE_LIST="$RNAQC_DIR/bam_list.txt"
find "$DATA_DIR" -name "*.bam" | while read -r bam; do
    sid=$(basename "$bam" | cut -d. -f1 | cut -d_ -f1 | tr '_' '-')
    echo "$bam $TARGET_BED $sid" >> "$SAMPLE_LIST"
done

# Run using xargs to manage 16 cores precisely
echo -e "externalsampleid\trna_skew" > "$SKEW_DATA"
cat "$SAMPLE_LIST" | xargs -P 16 -n 3 python3 "$RNAQC_DIR/worker.py" >> "$SKEW_DATA"

# ── STEP 2: METADATA, TISSUE, AND FINAL COVARIATES ────────────────────────────
echo "Merging results and generating final output tables..."

python3 - <<EOF
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter

# Remapping Logic
T_MAP = {"Motor_Cortex_Lateral": "Motor_Cortex", "Motor_Cortex_Medial": "Motor_Cortex", 
         "Lumbar_spinal_cord": "Lumbar_Spinal_Cord", "Cervical_spinal_cord": "Cervical_Spinal_Cord"}

# 1. Load Data
df = pd.read_csv("$METADATA")
df.columns = df.columns.str.strip()
wgs_meta = pd.read_csv("$WGS_META")
wgs_meta.columns = wgs_meta.columns.str.strip()

# 2. Tissue Tracking (Restored)
shared = sorted(set(wgs_meta["Externalsubjectid"]) & set(df["externalsubjectid"]))
patient_tissues = defaultdict(set)
sample_map = {row['externalsampleid']: row['externalsubjectid'] for _, row in df.iterrows()}
for d in Path("$DATA_DIR").glob("*/RNAseq/Processed/*/"):
    sid = d.name.replace("_", "-")
    if sid in sample_map and sample_map[sid] in shared:
        patient_tissues[sample_map[sid]].add(d.parts[-4])

with open("$PATIENT_TISSUES", "w") as f:
    f.write("subject_id\tn_tissues\ttissues\n")
    for s in sorted(patient_tissues.keys()):
        f.write(f"{s}\t{len(patient_tissues[s])}\t{'; '.join(sorted(patient_tissues[s]))}\n")

# 3. Covariate Cleaning
df['tissue'] = df['tissue'].map(lambda x: T_MAP.get(str(x).strip(), str(x).replace(" ", "_")))
df['sex'] = df['sex'].astype(str).str.capitalize()
mask = df['sex'].isin(['Unknown', 'Nd', 'Nan', ''])
df.loc[mask, 'sex'] = df[mask].apply(lambda r: 'Female' if 'XX' in str(r.get('sex_genotype','')) else 'Male', axis=1)

# Skew Rescue & Filtering
if Path("$SKEW_DATA").exists():
    sk = pd.read_csv("$SKEW_DATA", sep="\t").drop_duplicates('externalsampleid')
    df = df.merge(sk, on='externalsampleid', how='left')

df['rin_num'] = pd.to_numeric(df['rin'], errors='coerce')
df['sk_num'] = pd.to_numeric(df['rna_skew'], errors='coerce')
df['keep_for_qtl'] = df.apply(lambda r: True if pd.notna(r['rin_num']) else (-0.5 <= r['sk_num'] <= 0.5), axis=1)

# 4. Final Export
FULL_COLS = ["externalsampleid", "externalsubjectid", "sex", "subject_group", "age_at_death", "tissue", 
             "rin", "rna_skew", "keep_for_qtl", "disease_duration_in_months", "post_mortem_interval_in_hours", 
             "site_of_motor_onset", "c9orf72_repeat_expansion", "atxn2_repeat_expansion",
             "pct_african", "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas"]

df[[c for c in FULL_COLS if c in df.columns]].to_csv("$FINAL_COVARIATES", sep="\t", index=False)
EOF

echo "Process Complete."
