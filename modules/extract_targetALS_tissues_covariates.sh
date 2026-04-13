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
echo "Generating detailed summaries with cleaning and remapping..."

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

# Define Ancestry Columns here so they are available for cleaning and export
anc_cols = ["pct_african", "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas"]

# ── STEP 2: METADATA, TISSUE, AND FINAL COVARIATES ────────────────────────────
echo "Generating detailed summaries with cleaning and remapping..."

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

# Define Ancestry Columns
anc_cols = ["pct_african", "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas"]

# 2. REMAPPING DICTIONARIES
TISSUE_REMAP = {
    'Motor Cortex Lateral': 'Motor_Cortex', 'Motor Cortex Medial': 'Motor_Cortex',
    'Lateral Motor Cortex': 'Motor_Cortex', 'Medial Motor Cortex': 'Motor_Cortex',
    'Primary Motor Cortex L': 'Motor_Cortex', 'Primary Motor Cortex L ': 'Motor_Cortex',
    'Primary Motor Cortex M': 'Motor_Cortex', 'Cortex_Motor_BA4': 'Motor_Cortex', 
    'BA4 Motor Cortex': 'Motor_Cortex', 'Lateral motor cortex': 'Motor_Cortex', 
    'Cortex_Motor_Unspecified': 'Motor_Cortex', 'Motor Cortex': 'Motor_Cortex', 
    'Motor_Cortex': 'Motor_Cortex', 'Frontal Cortex': 'Frontal_Cortex', 
    'Cerebellum': 'Cerebellum', 'Cervical Spinal Cord': 'Cervical_Spinal_Cord', 
    'Cervical spinal cord': 'Cervical_Spinal_Cord', 'Spinal Cord Cervical': 'Cervical_Spinal_Cord', 
    'Spinal cord Cervical': 'Cervical_Spinal_Cord', 'Lumbar Spinal Cord': 'Lumbar_Spinal_Cord', 
    'Lumbar spinal cord': 'Lumbar_Spinal_Cord', 'Spinal Cord Lumbar': 'Lumbar_Spinal_Cord', 
    'Lumbosacral Spinal Cord': 'Lumbar_Spinal_Cord', 'Lumbosacaral spinal cord': 'Lumbar_Spinal_Cord', 
    'Spinal_Cord_Lumbosacral': 'Lumbar_Spinal_Cord', 'Thoracic Spinal Cord': 'Thoracic_Spinal_Cord'
}

SUBJECT_GROUP_REMAP = {
    'ALS Spectrum MND, Other Neurological Diseases': 'ALS Spectrum MND, Other Neurological Disorders',
    'Non Neurological Control': 'Non-Neurological Control',
    'Non-Neurological Control': 'Non-Neurological Control'
}

ONSET_REMAP = {
    'limb': 'Limb', 'Lower Limb': 'Limb', 'Upper Limb': 'Limb', 
    'Lower Extremity': 'Limb', 'Upper Extremity': 'Limb',
    'Lower limb': 'Limb', 'Upper limb': 'Limb', 'Upper and Lower Limbs': 'Limb',
    'Bulbar/Limb': 'Bulbar and Limb', 'Not Applicable': 'Not Applicable/NaN',
    'nan': 'Not Applicable/NaN'
}

# 3. CLEANING LOGIC
df = df.map(lambda x: x.strip() if isinstance(x, str) else x)

df['sex'] = df['sex'].replace({'male': 'Male', 'female': 'Female', 'ND': 'Unknown'})
df['subject_group'] = df['subject_group'].map(lambda x: SUBJECT_GROUP_REMAP.get(str(x), str(x)))
df['tissue'] = df['tissue'].map(lambda x: TISSUE_REMAP.get(str(x), str(x).replace(' ', '_')))
df['site_of_motor_onset'] = df['site_of_motor_onset'].map(lambda x: ONSET_REMAP.get(str(x), str(x)))

df['c9orf72_repeat_expansion'] = df['c9orf72_repeat_expansion'].replace({
    'yes': 'Yes', 'Negative': 'No', 'ND': 'Unknown', 'Not Applicable': 'Not Applicable/NaN'
})
df['c9orf72_repeat_expansion'] = df['c9orf72_repeat_expansion'].fillna('Not Applicable/NaN')

# Ancestry Cleaning and Binary Encoding
anc_binary_cols = []
for col in anc_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        df.loc[df[col] < 0, col] = np.nan
        
        # Create Binary Categorical Column (>= 50%)
        binary_name = "is_" + col.replace("pct_", "")
        df[binary_name] = (df[col] >= 0.50).astype(int)
        anc_binary_cols.append(binary_name)

# 4. QC Logic
df['rin_num'] = pd.to_numeric(df['rin'], errors='coerce')
df['sk_num'] = pd.to_numeric(df['rna_skew'], errors='coerce')
df['keep_for_qtl'] = df.apply(lambda r: True if (pd.notna(r['rin_num']) and r['rin_num'] > 5) or (-0.5 <= r['sk_num'] <= 0.5) else False, axis=1)

# Export Final Table
base_cols = ["externalsampleid", "externalsubjectid", "sex", "subject_group", "age_at_death", "tissue", 
             "rin", "rna_skew", "keep_for_qtl", "disease_duration_in_months", "post_mortem_interval_in_hours", 
             "site_of_motor_onset", "c9orf72_repeat_expansion"]
FULL_COLS = base_cols + anc_cols + anc_binary_cols

final_df = df[[c for c in FULL_COLS if c in df.columns]]
final_df.to_csv("$FINAL_COVARIATES", sep="\t", index=False)

# 5. GENERATE COVARIATE SUMMARY (Standard output continues as before)
with open("$COVARIATE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(f" target_ALS Covariate Summary\n Total samples: {len(df)} | Total subjects: {df['externalsubjectid'].nunique()}\n")
    f.write("============================================================\n\n")

    for col in ['sex', 'subject_group', 'tissue', 'site_of_motor_onset', 'c9orf72_repeat_expansion'] + anc_binary_cols:
        f.write(f"── {col.upper()} ──────────────────────────────\n")
        counts = df[col].value_counts(dropna=False)
        for val, count in counts.items():
            pct = (count / len(df)) * 100
            f.write(f"    {str(val):<45} {count:>5} ({pct:>4.1f}%)\n")
        f.write("\n")

# 6. TISSUE TRACKING (WGS + RNAseq)
wgs_subs = set(wgs_meta["Externalsubjectid"].dropna())
rna_subs = set(df["externalsubjectid"].dropna())
shared_subs = sorted(wgs_subs & rna_subs)

# Filter for shared subjects only
shared_df = df[df['externalsubjectid'].isin(shared_subs)]

# Subject mapping: {subject: [tissues]}
pt_tissues = shared_df.groupby('externalsubjectid')['tissue'].apply(lambda x: sorted(set(x))).to_dict()

with open("$PATIENT_TISSUES", "w") as f:
    f.write("subject_id\tn_tissues\ttissues\n")
    for sub in sorted(pt_tissues.keys()):
        ts = pt_tissues[sub]
        f.write(f"{sub}\t{len(ts)}\t{'; '.join(ts)}\n")

# TISSUE SUMMARY: Split by Subjects and Samples
# Calculate sample counts per tissue
sample_counts = shared_df['tissue'].value_counts().to_dict()
# Calculate unique subject counts per tissue
subject_counts = shared_df.groupby('tissue')['externalsubjectid'].nunique().to_dict()

with open("$TISSUE_SUMMARY", "w") as f:
    f.write(f"============================================================\n")
    f.write(f" target_ALS Tissue Summary (Shared WGS/RNA Patients)\n Shared Patients: {len(shared_subs)}\n")
    f.write(f"============================================================\n\n")
    f.write(f"{'TISSUE':<35} {'SUBJECTS':<12} {'SAMPLES':<10}\n")
    f.write(f"{'-'*34:<35} {'-'*10:<12} {'-'*8:<10}\n")
    
    # Sort by subject count descending
    for t in sorted(subject_counts.keys(), key=lambda x: -subject_counts[x]):
        f.write(f"{str(t):<35} {subject_counts[t]:<12} {sample_counts[t]:<10}\n")
EOF
