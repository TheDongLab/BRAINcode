#!/bin/bash
#SBATCH --job-name=targetALS_complete_audit
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --time=23:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%j.err

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

# Helper for robust remapping
def clean_value(val, mapping):
    if pd.isna(val): return "Not Applicable/NaN"
    clean = str(val).strip().lower().replace(" ", "_")
    return mapping.get(clean, str(val).strip())

# 2. DEFINITIVE REMAPPING DICTIONARIES
T_REMAP = {
    'motor_cortex_lateral': 'Motor_Cortex', 'motor_cortex_medial': 'Motor_Cortex',
    'lateral_motor_cortex': 'Motor_Cortex', 'medial_motor_cortex': 'Motor_Cortex',
    'primary_motor_cortex_l': 'Motor_Cortex', 'primary_motor_cortex_m': 'Motor_Cortex',
    'cortex_motor_ba4': 'Motor_Cortex', 'ba4_motor_cortex': 'Motor_Cortex',
    'cortex_motor_unspecified': 'Motor_Cortex', 'motor_cortex': 'Motor_Cortex',
    'frontal_cortex': 'Frontal_Cortex', 'cerebellum': 'Cerebellum',
    'cervical_spinal_cord': 'Cervical_Spinal_Cord', 'spinal_cord_cervical': 'Cervical_Spinal_Cord',
    'lumbar_spinal_cord': 'Lumbar_Spinal_Cord', 'spinal_cord_lumbar': 'Lumbar_Spinal_Cord',
    'lumbosacral_spinal_cord': 'Lumbar_Spinal_Cord', 'lumbosacaral_spinal_cord': 'Lumbar_Spinal_Cord',
    'thoracic_spinal_cord': 'Thoracic_Spinal_Cord'
}

G_REMAP = {
    'other_mnd': 'Other Neurological Disorders',
    'other_neurological_disorders': 'Other Neurological Disorders',
    'als_spectrum_mnd,_other_neurological_diseases': 'ALS Spectrum MND, Other Neurological Disorders',
    'als_spectrum_mnd,_other_neurological_disorders': 'ALS Spectrum MND, Other Neurological Disorders',
    'alzheimer’s_disease,_definite:_cerad_criteria,ftld-mnd/mni,amyotrophic_lateral_sclerosis': 'ALS/FTLD Spectrum',
    'ftld-mnd/mni,amyotrophic_lateral_sclerosis': 'ALS/FTLD Spectrum',
    'ftld-tdp,_cerebrovascular_disease': 'ALS/FTLD Spectrum',
    'ftd,_tdp43_subtype': 'ALS/FTLD Spectrum',
    'ftld': 'ALS/FTLD Spectrum',
    'ftld_with_tau_positive_inclusions': 'ALS/FTLD Spectrum',
    'non_neurological_control': 'Non-Neurological Control'
}

O_REMAP = {
    'limb': 'Limb', 'lower_limb': 'Limb', 'upper_limb': 'Limb', 
    'lower_extremity': 'Limb', 'upper_extremity': 'Limb',
    'upper_and_lower_limbs': 'Limb', 'bulbar/limb': 'Bulbar and Limb',
    'not_applicable': 'Not Applicable/NaN'
}

C_REMAP = {
    'yes': 'Yes', 'negative': 'No', 'nd': 'Unknown', 
    'not_applicable': 'Not Applicable/NaN', 'abnormal_repeat_expansion': 'Yes'
}

# 3. GLOBAL CLEANING
df['sex'] = df['sex'].str.strip().str.capitalize().replace({'Male': 'Male', 'Female': 'Female', 'Nd': 'Unknown'})
df['tissue'] = df['tissue'].apply(lambda x: clean_value(x, T_REMAP))
df['subject_group'] = df['subject_group'].apply(lambda x: clean_value(x, G_REMAP))
df['site_of_motor_onset'] = df['site_of_motor_onset'].apply(lambda x: clean_value(x, O_REMAP))
df['c9orf72_repeat_expansion'] = df['c9orf72_repeat_expansion'].apply(lambda x: clean_value(x, C_REMAP))

# ── ANCESTRY FIX ──────────────────────────────────────────────────────────────
anc_cols = ["pct_african", "pct_south_asian", "pct_east_asian", "pct_european", "pct_americas"]
anc_binary_cols = []
for col in anc_cols:
    if col in df.columns:
        # Force numeric conversion in case strip/map turned them into objects
        df[col] = pd.to_numeric(df[col], errors='coerce')
        bin_name = "is_" + col.replace("pct_", "")
        # Create binary flag: 1 if >= 0.5, 0 otherwise (NaNs become 0)
        df[bin_name] = (df[col] >= 0.5).astype(int)
        anc_binary_cols.append(bin_name)
# ──────────────────────────────────────────────────────────────────────────────

# Split Cohorts
wgs_subs = set(wgs_meta["Externalsubjectid"].dropna())
shared_df = df[df['externalsubjectid'].isin(wgs_subs)]

# 4. EXPORT MAIN TABLE
df.to_csv("$FINAL_COVARIATES", sep="\t", index=False)

# 5. COVARIATE SUMMARY (Audit)
with open("$COVARIATE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(" target_ALS Covariate Frequency Audit\n")
    f.write(f" GLOBAL (All RNA):             {len(df):<5} samples | {df['externalsubjectid'].nunique():<5} subjects\n")
    f.write(f" SHARED (also has WGS data):   {len(shared_df):<5} samples | {shared_df['externalsubjectid'].nunique():<5} subjects\n")
    f.write("============================================================\n\n")

    cat_cols = ['sex', 'subject_group', 'tissue', 'site_of_motor_onset', 'c9orf72_repeat_expansion'] + anc_binary_cols
    for col in cat_cols:
        if col not in df.columns: continue
        f.write(f"── {col.upper():<35} {'GLOBAL':<15} {'SHARED (also has WGS data)':<15}\n")
        g_counts = df[col].value_counts(dropna=False)
        s_counts = shared_df[col].value_counts(dropna=False)
        all_labels = sorted(set(g_counts.index.astype(str)) | set(s_counts.index.astype(str)))
        for val in all_labels:
            g_c, s_c = g_counts.get(val, 0), s_counts.get(val, 0)
            f.write(f"    {val:<31} {g_c:>5} ({ (g_c/len(df))*100:>4.1f}%)    {s_c:>5} ({ (s_c/len(shared_df))*100 if len(shared_df)>0 else 0:>4.1f}%)\n")
        f.write("\n")

    f.write("── ANCESTRY PERCENTAGE BINS ───────────────────────────\n\n")
    bins = [0, 0.01, 0.05, 0.25, 0.50, 0.75, 1.01]
    labels = ["<1%", "1–5%", "5–25%", "25–50%", "50–75%", "75–100%"]
    for col in anc_cols:
        if col not in df.columns: continue
        f.write(f"── {col.upper():<35} {'GLOBAL':<15} {'SHARED (also has WGS data)':<15}\n")
        g_valid = df[col].dropna()
        s_valid = shared_df[col].dropna()
        g_binned = pd.cut(g_valid, bins=bins, labels=labels, right=False).value_counts().sort_index()
        s_binned = pd.cut(s_valid, bins=bins, labels=labels, right=False).value_counts().sort_index()
        for b in labels:
            gc, sc = g_binned.get(b, 0), s_binned.get(b, 0)
            f.write(f"    {b:<31} {gc:>5} ({ (gc/len(df))*100:>4.1f}%)    {sc:>5} ({ (sc/len(shared_df))*100 if len(shared_df)>0 else 0:>4.1f}%)\n")
        f.write("\n")

# 6. TISSUE SUMMARY
subject_counts = shared_df.groupby('tissue')['externalsubjectid'].nunique()
sample_counts = shared_df['tissue'].value_counts()
with open("$TISSUE_SUMMARY", "w") as f:
    f.write(f"============================================================\n")
    f.write(f" target_ALS Tissue Summary (Shared Patients Only)\n")
    f.write(f" Shared Patients Total: {shared_df['externalsubjectid'].nunique()}\n")
    f.write(f"============================================================\n\n")
    f.write(f"{'TISSUE':<35} {'SUBJECTS':<15} {'SAMPLES':<15}\n")
    f.write(f"{'-'*34:<35} {'-'*12:<15} {'-'*11:<15}\n")
    for t in sorted(subject_counts.index, key=lambda x: -subject_counts[x]):
        f.write(f"{str(t):<35} {subject_counts[t]:<15} {sample_counts[t]:<15}\n")

# 7. PATIENT TISSUE BREAKDOWN
pt_tissues = shared_df.groupby('externalsubjectid')['tissue'].apply(lambda x: sorted(set(x))).to_dict()
with open("$PATIENT_TISSUES", "w") as f:
    f.write("subject_id\tn_tissues\ttissues\n")
    for sub in sorted(pt_tissues.keys()):
        ts = pt_tissues[sub]
        f.write(f"{sub}\t{len(ts)}\t{'; '.join(ts)}\n")
EOF

echo "Done. Complete audit outputs in $OUTDIR"
