#!/bin/bash
#SBATCH --job-name=targetALS_covariate_tissue_summary
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=23:00:00
#SBATCH -p day

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
if [ ! -f "$TARGET_BED" ]; then
    echo 'Cleaning BED file and selecting longest genes...'
    set +o pipefail
    grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)[[:space:]]' "$ORIG_BED" | \
    awk -v OFS='\t' '{split($4, a, "___"); $4=a[1]; print $1, $2, $3, $4, $5, $6, $3-$2}' | \
    sort -k7,7rn | \
    head -n 500 | \
    cut -f1-6 > "$TARGET_BED" || true
    set -o pipefail
fi

# ── STEP 1: PARALLEL SKEWNESS CALCULATION (WITH SKIP LOGIC) ───────────────────
if [ ! -f "$SKEW_DATA" ]; then
    echo -e "externalsampleid\trna_skew" > "$SKEW_DATA"
fi

PROCESSED_IDS=$(tail -n +2 "$SKEW_DATA" | cut -f1)
> "$SAMPLE_LIST"

find "$DATA_DIR" -name "STAR.Aligned.sortedByCoord.out.bam" | while read -r bam; do
    sid=$(basename "$(dirname "$bam")" | tr '_' '-')
    if ! echo "$PROCESSED_IDS" | grep -qxw "$sid"; then
        echo "$bam $TARGET_BED $sid" >> "$SAMPLE_LIST"
    fi
done

N_TO_DO=$(wc -l < "$SAMPLE_LIST")
echo "Already processed: $(echo "$PROCESSED_IDS" | wc -l). Remaining: $N_TO_DO"

if [ "$N_TO_DO" -gt 0 ]; then
    cat << 'EOF' > "$RNAQC_DIR/worker.py"
import sys, subprocess, numpy as np
def compute(bam, bed):
    gene_lookup = {}
    with open(bed) as f:
        for line in f:
            p = line.split()
            chrom, start, end = p[0], int(p[1]), int(p[2])
            if chrom not in gene_lookup: gene_lookup[chrom] = []
            gene_lookup[chrom].append((start, end))
    cmd = ["samtools", "depth", "-a", "-b", bed, bam]
    try:
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc:
            bins = np.zeros(101)
            for line in proc.stdout:
                parts = line.split()
                if len(parts) < 3: continue
                chrom, pos, depth = parts[0], int(parts[1]), int(parts[2])
                if chrom in gene_lookup:
                    for g_start, g_end in gene_lookup[chrom]:
                        if g_start <= pos <= g_end:
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
    print(f"{sys.argv[3]}\t{compute(sys.argv[1], sys.argv[2])}", flush=True)
EOF
    cat "$SAMPLE_LIST" | xargs -P 16 -n 3 python3 "$RNAQC_DIR/worker.py" >> "$SKEW_DATA"
fi

# ── STEP 2: METADATA, TISSUE, AND FINAL COVARIATES ────────────────────────────
echo "Generating detailed summaries..."

python3 - <<EOF
import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter

# 1. Load Data
df = pd.read_csv("$METADATA")
wgs_meta = pd.read_csv("$WGS_META")
sk = pd.read_csv("$SKEW_DATA", sep="\t").drop_duplicates('externalsampleid')
df = df.merge(sk, on='externalsampleid', how='left')

# 2. Cleanup & QC Logic
df['rin_num'] = pd.to_numeric(df['rin'], errors='coerce')
df['sk_num'] = pd.to_numeric(df['rna_skew'], errors='coerce')
df['keep_for_qtl'] = df.apply(lambda r: True if (pd.notna(r['rin_num']) and r['rin_num'] > 5) or (-0.5 <= r['sk_num'] <= 0.5) else False, axis=1)

FULL_COLS = ["externalsampleid", "externalsubjectid", "sex", "subject_group", "age_at_death", "tissue", 
             "rin", "rna_skew", "keep_for_qtl", "disease_duration_in_months", "post_mortem_interval_in_hours", 
             "site_of_motor_onset", "c9orf72_repeat_expansion", "pct_african", "pct_south_asian", 
             "pct_east_asian", "pct_european", "pct_americas"]

final_df = df[[c for c in FULL_COLS if c in df.columns]]
final_df.to_csv("$FINAL_COVARIATES", sep="\t", index=False)

# 3. Covariate Summary Report
with open("$COVARIATE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(f" target_ALS Covariate Summary\n Total samples: {len(df)} | Total subjects: {df['externalsubjectid'].nunique()}\n")
    f.write("============================================================\n\n")

    for col in ['sex', 'subject_group', 'tissue', 'site_of_motor_onset', 'c9orf72_repeat_expansion']:
        f.write(f"── {col.upper()} ──────────────────────────────\n")
        counts = df[col].value_counts(dropna=False)
        for val, count in counts.items():
            pct = (count / len(df)) * 100
            f.write(f"    {str(val):<45} {count:>5} ({pct:>4.1f}%)\n")
        f.write("\n")

    f.write("── NUMERICAL COVARIATES ────────────────────────────────\n\n")
    for col in ['age_at_death', 'rin', 'disease_duration_in_months', 'post_mortem_interval_in_hours']:
        data = pd.to_numeric(df[col], errors='coerce').dropna()
        f.write(f"{col} (n={len(data)}, {len(df)-len(data)} missing)\n")
        f.write(f"    mean: {data.mean():.2f} | median: {data.median():.2f} | std: {data.std():.2f}\n\n")

    f.write("── ANCESTRY ───────────────────────────────────────────\n\n")
    anc_cols = ["pct_african", "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas"]
    bins = [0, 0.01, 0.05, 0.25, 0.50, 0.75, 1.01]
    labels = ["<1%", "1–5%", "5–25%", "25–50%", "50–75%", "75–100%"]
    
    for col in anc_cols:
        data = pd.to_numeric(df[col], errors='coerce').dropna()
        f.write(f"{col} (n={len(data)}, {len(df)-len(data)} missing/invalid)\n")
        f.write(f"    mean: {data.mean():.4f} | median: {data.median():.4f}\n")
        binned = pd.cut(data, bins=bins, labels=labels, right=False).value_counts().sort_index()
        for b, c in binned.items():
            pct = (c / len(data)) * 100 if len(data)>0 else 0
            f.write(f"      {b:<15} {c:>5} ({pct:>4.1f}%)\n")
        f.write("\n")

# 4. Tissue Tracking & Shared Summary
wgs_subs = set(wgs_meta["Externalsubjectid"].dropna())
rna_subs = set(df["externalsubjectid"].dropna())
shared = sorted(wgs_subs & rna_subs)

# Map subjects to tissues for those in shared set
shared_df = df[df['externalsubjectid'].isin(shared)]
pt_tissues = shared_df.groupby('externalsubjectid')['tissue'].apply(lambda x: sorted(set(x))).to_dict()

with open("$PATIENT_TISSUES", "w") as f:
    f.write("subject_id\tn_tissues\ttissues\n")
    for sub in sorted(pt_tissues.keys()):
        ts = pt_tissues[sub]
        f.write(f"{sub}\t{len(ts)}\t{'; '.join(ts)}\n")

t_counts = Counter([t for ts in pt_tissues.values() for t in ts])
with open("$TISSUE_SUMMARY", "w") as f:
    f.write(f"==============================\n target_ALS Tissue Summary\n Shared Patients: {len(shared)}\n==============================\n\n")
    f.write(f"{'COUNT':<8} {'TISSUE'}\n")
    f.write(f"{'-----':<8} {'------------------------------'}\n")
    for t, c in sorted(t_counts.items(), key=lambda x: -x[1]):
        f.write(f"{c:<8} {str(t).replace(' ', '_')}\n")
EOF

echo "Process Complete."
echo "Summary available at: $COVARIATE_SUMMARY"
