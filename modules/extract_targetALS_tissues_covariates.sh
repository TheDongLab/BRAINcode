#!/bin/bash
#SBATCH --job-name=targetALS_covariate_tissue_summary
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=20:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/targetALS_covariate_tissue_summary.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/targetALS_covariate_tissue_summary.err

# ── target_ALS tissue summary + covariate extraction ─────────────────────────
# Optimized to target STAR symlinks and handle global tissues
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
export TARGET_BED="$RNAQC_DIR/reference_subset.bed"
SKEW_DATA="$RNAQC_DIR/calculated_skew.tsv"
FINAL_COVARIATES="$OUTDIR/covariates.tsv"
PATIENT_TISSUES="$OUTDIR/patient_tissue_breakdown.tsv"
TISSUE_SUMMARY="$OUTDIR/tissue_summary.txt"
COVARIATE_SUMMARY="$OUTDIR/covariate_summary.txt"
SAMPLE_LIST="$RNAQC_DIR/bam_list.txt"

# ── STEP 0: BED REFORMATTING ──────────────────────────────────────────────────
echo 'Cleaning BED file and selecting longest genes...'

# Prevent SIGPIPE from head -n 500 killing the script
set +o pipefail
grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]' "$ORIG_BED" | \
awk -v OFS='\t' '{split($4, a, "___"); $4=a[1]; print $1, $2, $3, $4, $5, $6, $3-$2}' | \
sort -k7,7rn | \
head -n 500 | \
cut -f1-6 > "$TARGET_BED" || true
set -o pipefail

# ── STEP 1: PARALLEL SKEWNESS CALCULATION ─────────────────────────────────────
echo "Generating filtered BAM list targeting STAR symlinks..."

# Clear list and find only the specific STAR symlinks across all tissue folders
> "$SAMPLE_LIST"
find "$DATA_DIR" -name "STAR.Aligned.sortedByCoord.out.bam" | while read -r bam; do
    # Get folder name (e.g. SD_001_23_MCX) and flip underscores to dashes for metadata matching
    sid=$(basename "$(dirname "$bam")" | tr '_' '-')
    echo "$bam $TARGET_BED $sid" >> "$SAMPLE_LIST"
done

echo "Found $(wc -l < "$SAMPLE_LIST") STAR BAMs. Starting Skew Calculation..."

# Worker script - optimized for memory and speed
cat << 'EOF' > "$RNAQC_DIR/worker.py"
import sys, subprocess, numpy as np
def compute(bam, bed):
    cmd = ["samtools", "depth", "-a", "-b", bed, bam]
    try:
        # Stream the depth data to avoid memory bloat
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc:
            bins = np.zeros(101)
            genes = []
            with open(bed) as f:
                for line in f:
                    p = line.split()
                    genes.append((p[0], int(p[1]), int(p[2])))
            
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
    except Exception: return "NaN"

if __name__ == "__main__":
    bam_path, bed_path, sid = sys.argv[1], sys.argv[2], sys.argv[3]
    print(f"{sid}\t{compute(bam_path, bed_path)}")
EOF

# Run parallel using xargs
echo -e "externalsampleid\trna_skew" > "$SKEW_DATA"
cat "$SAMPLE_LIST" | xargs -P 16 -n 3 python3 "$RNAQC_DIR/worker.py" >> "$SKEW_DATA"

# ── STEP 2: METADATA, TISSUE, AND FINAL COVARIATES ────────────────────────────
echo "Merging skew results with $METADATA..."

python3 - <<EOF
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

# 1. Load Data
df = pd.read_csv("$METADATA")
df.columns = df.columns.str.strip()
wgs_meta = pd.read_csv("$WGS_META")
wgs_meta.columns = wgs_meta.columns.str.strip()

# 2. Skew Integration
if Path("$SKEW_DATA").exists():
    sk = pd.read_csv("$SKEW_DATA", sep="\t").drop_duplicates('externalsampleid')
    # Ensure sample IDs match by stripping potential trailing info
    df = df.merge(sk, on='externalsampleid', how='left')

# 3. Covariate Cleaning & Rescue Logic
# If RIN is missing, we use RNA Skew as a proxy for degradation
df['rin_num'] = pd.to_numeric(df['rin'], errors='coerce')
df['sk_num'] = pd.to_numeric(df['rna_skew'], errors='coerce')

# Criteria: Keep if RIN exists OR if Skew is within healthy bounds (-0.5 to 0.5)
df['keep_for_qtl'] = df.apply(
    lambda r: True if pd.notna(r['rin_num']) and r['rin_num'] > 5 
    else (True if (-0.5 <= r['sk_num'] <= 0.5) else False), 
    axis=1
)

# 4. Final Table Formatting
FULL_COLS = ["externalsampleid", "externalsubjectid", "sex", "subject_group", "age_at_death", "tissue", 
             "rin", "rna_skew", "keep_for_qtl", "disease_duration_in_months", "post_mortem_interval_in_hours", 
             "site_of_motor_onset", "c9orf72_repeat_expansion", "atxn2_repeat_expansion",
             "pct_african", "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas"]

# Filter to available columns and export
final_df = df[[c for c in FULL_COLS if c in df.columns]]
final_df.to_csv("$FINAL_COVARIATES", sep="\t", index=False)

# 5. Tissue Breakdown Report
shared_subjects = set(wgs_meta["Externalsubjectid"]) & set(df["externalsubjectid"])
summary = df[df['externalsubjectid'].isin(shared_subjects)].groupby('externalsubjectid')['tissue'].apply(lambda x: '; '.join(sorted(set(x)))).reset_index()
summary['n_tissues'] = summary['tissue'].str.count(';') + 1
summary.columns = ['subject_id', 'tissues', 'n_tissues']
summary[['subject_id', 'n_tissues', 'tissues']].to_csv("$PATIENT_TISSUES", sep="\t", index=False)
EOF

echo "Process Complete. Final covariates: $FINAL_COVARIATES"
