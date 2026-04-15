#!/bin/bash
#SBATCH --job-name=extract_targetALS_tissues_covariates
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --time=23:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/extract_targetALS_tissues_covariates.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/extract_targetALS_tissues_covariates.err

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

# ── STEP 1: PARALLEL SKEWNESS CALCULATION ─────────────────────────────────────
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

# ── STEP 2: METADATA CONSOLIDATION AND AUDIT ──────────────────────────────────
echo "Consolidating metadata and generating complete audit..."

python3 - <<EOF
import pandas as pd
import numpy as np

# 1. LOAD DATA
df = pd.read_csv("$METADATA")
wgs_meta = pd.read_csv("$WGS_META")
sk = pd.read_csv("$SKEW_DATA", sep="\t").drop_duplicates('externalsampleid')
df = df.merge(sk, on='externalsampleid', how='left')

anc_cols = ["pct_african", "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas"]

# 2. DEFINITIVE REMAPPING
TISSUE_REMAP = {
    'Motor Cortex Lateral': 'Motor_Cortex', 'Motor Cortex Medial': 'Motor_Cortex',
    'Lateral Motor Cortex': 'Motor_Cortex', 'Medial Motor Cortex': 'Motor_Cortex',
    'Primary Motor Cortex L': 'Motor_Cortex', 'Primary Motor Cortex M': 'Motor_Cortex', 'Cortex_Motor_Unspecified': 'Motor_Cortex',
    'Cortex_Motor_BA4': 'Motor_Cortex', 'BA4 Motor Cortex': 'Motor_Cortex', 'Lateral_motor_cortex': 'Motor_Cortex',
    'Frontal Cortex': 'Frontal_Cortex', 'Cerebellum': 'Cerebellum', 'Spinal_Cord_Cervical': 'Cervical_Spinal_Cord', 
    'Cervical Spinal Cord': 'Cervical_Spinal_Cord', 'Cervical_spinal_cord': 'Cervical_Spinal_Cord', 'Spinal_cord_Cervical': 'Cervical_Spinal_Cord', 'Lumbar Spinal Cord': 'Lumbar_Spinal_Cord',
    'Thoracic Spinal Cord': 'Thoracic_Spinal_Cord', 'Spinal_Cord_Lumbosacral': 'Lumbar_Spinal_Cord', 'Lumbosacral_Spinal_Cord': 'Lumbar_Spinal_Cord', 'Lumbar_spinal_cord': 'Lumbar_Spinal_Cord'
}

SUBJECT_GROUP_REMAP = {
    'Other MND': 'Other Neurological Disorders',
    'Other Neurological Disorders': 'Other Neurological Disorders',
    'ALS Spectrum MND, Other Neurological Diseases': 'ALS Spectrum MND, Other Neurological Disorders',
    'Alzheimer’s Disease, Definite: CERAD criteria,FTLD-MND/MNI,Amyotrophic Lateral Sclerosis': 'ALS/FTLD Spectrum',
    'FTLD-MND/MNI,Amyotrophic Lateral Sclerosis': 'ALS/FTLD Spectrum',
    'FTLD-TDP, Cerebrovascular disease': 'ALS/FTLD Spectrum',
    'FTD, TDP43 subtype': 'ALS/FTLD Spectrum',
    'FTLD': 'ALS/FTLD Spectrum',
    'FTLD with tau positive inclusions': 'ALS/FTLD Spectrum',
    'Non Neurological Control': 'Non-Neurological Control'
}

ONSET_REMAP = {
    'limb': 'Limb', 'Lower Limb': 'Limb', 'Upper Limb': 'Limb', 
    'Lower Extremity': 'Limb', 'Upper Extremity': 'Limb',
    'Lower limb': 'Limb', 'Upper limb': 'Limb', 'Upper and Lower Limbs': 'Limb',
    'Bulbar/Limb': 'Bulbar and Limb', 'Not Applicable': 'Not Applicable/NaN',
    'nan': 'Not Applicable/NaN', 'Unknown': 'Not Applicable/NaN'
}

C9_REMAP = {
    'yes': 'Yes', 'Negative': 'No', 'ND': 'Unknown', 
    'Not Applicable': 'Not Applicable/NaN', 'Abnormal Repeat Expansion': 'Yes',
    'nan': 'Not Applicable/NaN', 'Unknown': 'Not Applicable/NaN'
}

# 3. GLOBAL CLEANING (Preserving logic exactly as requested)
df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
df['sex'] = df['sex'].replace({'male': 'Male', 'female': 'Female', 'ND': 'Unknown'})
df['tissue'] = df['tissue'].map(lambda x: TISSUE_REMAP.get(str(x), str(x).replace(' ', '_')))
df['subject_group'] = df['subject_group'].map(lambda x: SUBJECT_GROUP_REMAP.get(str(x), str(x)))
df['site_of_motor_onset'] = df['site_of_motor_onset'].map(lambda x: ONSET_REMAP.get(str(x), str(x)))
df['c9orf72_repeat_expansion'] = df['c9orf72_repeat_expansion'].map(lambda x: C9_REMAP.get(str(x), str(x)))
df['c9orf72_repeat_expansion'] = df['c9orf72_repeat_expansion'].fillna('Not Applicable/NaN')

anc_binary_cols = []
for col in anc_cols:
    if col in df.columns:
        # Re-establishing numeric type specifically for these columns
        df[col] = pd.to_numeric(df[col], errors='coerce')
        df.loc[df[col] < 0, col] = np.nan
        bin_name = "is_" + col.replace("pct_", "")
        df[bin_name] = (df[col] >= 0.50).astype(int)
        anc_binary_cols.append(bin_name)

# Split Cohorts
wgs_subs = set(wgs_meta["Externalsubjectid"].dropna())
shared_df = df[df['externalsubjectid'].isin(wgs_subs)].copy()
# FORCE tissue column normalization AGAIN
shared_df['tissue'] = shared_df['tissue'].astype(str)

# 4. EXPORT MAIN TABLE
df.to_csv("$FINAL_COVARIATES", sep="\t", index=False)

# 5. COVARIATE SUMMARY (Audit)
with open("$COVARIATE_SUMMARY", "w") as f:

    f.write("============================================================\n")
    f.write(" target_ALS Covariate Frequency Audit\n")
    f.write(f" GLOBAL (All RNA):             {len(df):<5} samples | {df['externalsubjectid'].nunique():<5} subjects\n")
    f.write(f" SHARED (also has WGS data):   {len(shared_df):<5} samples | {shared_df['externalsubjectid'].nunique():<5} subjects\n")
    f.write("============================================================\n\n")

    cat_cols = ['sex', 'subject_group', 'site_of_motor_onset', 'c9orf72_repeat_expansion'] + anc_binary_cols
    for col in cat_cols:
        if col not in df.columns: continue
        f.write(f"── {col.upper():<35} {'GLOBAL':<15} {'SHARED (also has WGS data)':<15}\n")
        g_counts = df[col].value_counts(dropna=False)
        s_counts = shared_df[col].value_counts(dropna=False)
        all_labels = sorted(set(g_counts.index.astype(str)) | set(s_counts.index.astype(str)))
        for val in all_labels:
            orig_val = next((k for k in g_counts.index if str(k) == val), val)
            g_c, s_c = g_counts.get(orig_val, 0), s_counts.get(orig_val, 0)
            f.write(f"    {val:<31} {g_c:>5} ({ (g_c/len(df))*100:>4.1f}%)    {s_c:>5} ({ (s_c/len(shared_df))*100 if len(shared_df)>0 else 0:>4.1f}%)\n")
        f.write("\n")

    f.write("── ANCESTRY PERCENTAGE BINS ───────────────────────────\n\n")
    bins = [0, 0.01, 0.05, 0.25, 0.50, 0.75, 1.01]
    labels = ["<1%", "1–5%", "5–25%", "25–50%", "50–75%", "75–100%"]
    for col in anc_cols:
        f.write(f"── {col.upper():<35} {'GLOBAL':<15} {'SHARED (also has WGS data)':<15}\n")
        g_data, s_data = df[col].dropna(), shared_df[col].dropna()
        g_binned = pd.cut(g_data, bins=bins, labels=labels, right=False).value_counts().sort_index()
        s_binned = pd.cut(s_data, bins=bins, labels=labels, right=False).value_counts().sort_index()
        for b in labels:
            gc, sc = g_binned.get(b, 0), s_binned.get(b, 0)
            f.write(f"    {b:<31} {gc:>5} ({ (gc/len(g_data))*100 if len(g_data)>0 else 0:>4.1f}%)    {sc:>5} ({ (sc/len(s_data))*100 if len(s_data)>0 else 0:>4.1f}%)\n")
        f.write("\n")

    f.write("── NUMERICAL COVARIATES ────────────────────────────────────\n\n")
    num_cols = {
        'age_at_death': 'age_at_death', 
        'rin': 'rin', 
        'disease_duration_in_months': 'disease_duration_in_months', 
        'post_mortem_interval_in_hours': 'post_mortem_interval_in_hours'
    }
    for label, col in num_cols.items():
        if col in df.columns:
            vals = pd.to_numeric(df[col], errors='coerce').dropna()
            missing = len(df) - len(vals)
            f.write(f"{label} (n={len(vals)}, {missing} missing)\n")
            f.write(f"    mean: {vals.mean():.2f} | median: {vals.median():.2f} | std: {vals.std():.2f}\n\n")
    # -------------------------------------

# 6. TISSUE SUMMARY
# enforce remap BEFORE any counting
shared_df['tissue'] = shared_df['tissue'].map(
    lambda x: TISSUE_REMAP.get(x, x)
)

subject_counts = shared_df.groupby('tissue')['externalsubjectid'].nunique()
sample_counts = shared_df['tissue'].value_counts()

with open("$TISSUE_SUMMARY", "w") as f:
    f.write(f"============================================================\n")
    f.write(f" target_ALS Tissue Summary (Shared Patients Only)\n")
    f.write(f" Shared Patients: {len(wgs_subs & set(df['externalsubjectid']))}\n")
    f.write(f"============================================================\n\n")
    f.write(f"{'TISSUE':<35} {'SUBJECTS':<15} {'SAMPLES':<15}\n")
    f.write(f"{'-'*34:<35} {'-'*12:<15} {'-'*11:<15}\n")
    for t in sorted(subject_counts.index, key=lambda x: -subject_counts[x]):
        f.write(f"{str(t):<35} {subject_counts[t]:<15} {sample_counts[t]:<15}" + "\n")


# OPTIONAL SAFETY (prevents silent mismatches like trailing spaces)
shared_df['tissue'] = shared_df['tissue'].str.strip()

# 7. PATIENT TISSUE BREAKDOWN
pt_tissues = shared_df.groupby('externalsubjectid')['tissue'].apply(lambda x: sorted(set(x))).to_dict()
with open("$PATIENT_TISSUES", "w") as f:
    f.write("subject_id\tn_tissues\ttissues\n")
    for sub in sorted(pt_tissues.keys()):
        ts = pt_tissues[sub]
        f.write(f"{sub}\t{len(ts)}\t{'; '.join(ts)}\n")
EOF

echo "Done. Complete audit outputs in $OUTDIR"
