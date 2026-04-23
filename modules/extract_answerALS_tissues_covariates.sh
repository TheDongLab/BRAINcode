#!/bin/bash
#SBATCH --job-name=extract_answerALS_covariates
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/extract_answerALS_audit.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/extract_answerALS_audit.err

set -euo pipefail

# ── PATHS ─────────────────────────────────────────────────────────────────────
# Executing from /home/zw529/donglab/data/answer_ALS/QTL
DATA_DIR="/home/zw529/donglab/data/answer_ALS"
OUTDIR="$DATA_DIR/QTL"
CONSOLIDATED_META="$DATA_DIR/consolidated_metadata_full.tsv"

# Ensure the output directory exists
mkdir -p "$OUTDIR"

# Output Files
COVARIATE_SUMMARY="$OUTDIR/covariate_summary.txt"
TISSUE_SUMMARY="$OUTDIR/tissue_summary.txt"
DIST_PLOT="$OUTDIR/answerALS_distributions.png"
CORR_PLOT="$OUTDIR/answerALS_collinearity_heatmap.png"

# ── STEP 1: SUBJECT-LEVEL AGGREGATION & AUDIT ─────────────────────────────────
echo "Collapsing metadata to Subject-level and generating audit..."

python3 - <<EOF
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# 1. LOAD DATA
df = pd.read_csv("$CONSOLIDATED_META", sep="\t", low_memory=False)

# 2. DEFINE AGGREGATION LOGIC
# For clinical/categorical info, we take the 'first' or 'max'
# For biological measurements (staining), we take the 'mean' per subject
agg_dict = {
    'SEX': 'first',
    'SUBJECT_GROUP': 'first',
    'SUBJECT_SUBGROUP': 'first',
    'Site_of_Onset': 'first',
    'EH_C9orf72': 'max',
    'has_discordant_sample': 'max',
    'AGE_AT_SYMPTOM_ONSET': 'max',
    'AGE_AT_DEATH': 'max',
    'Age_Sample_Collection_Cedars': 'max',
    'Age_at_First_PBMC_Collection': 'max',
    'ALSFRS_R_Baseline_Value': 'max',
    'ALSFRS_R_Latest_Value': 'min', # Getting the most recent/lowest score
    'ALSFRS_R_PROGRESSION_SLOPE': 'mean',
    'C9orf72_repeat_length': 'max',
    'ATXN2_repeat_length': 'max',
    'PCT_SMI32': 'mean',
    'PCT_ISL1': 'mean',
    'PCT_NKX61': 'mean',
    'PCT_TUJ1': 'mean',
    'PCT_S100b': 'mean',
    'PCT_Nestin': 'mean'
}

# Only aggregate columns that actually exist in the dataframe
agg_dict = {k: v for k, v in agg_dict.items() if k in df.columns}

# Perform the collapse
df_sub = df.groupby('Participant_ID').agg(agg_dict).reset_index()

# 3. GENERATE FIGURES (Subject-Level)
print(f"Generating figures for N={len(df_sub)} subjects...")
age_cols = ['AGE_AT_SYMPTOM_ONSET', 'AGE_AT_DEATH', 'Age_Sample_Collection_Cedars', 'Age_at_First_PBMC_Collection']
ipsc_qc = ['PCT_SMI32', 'PCT_ISL1', 'PCT_NKX61', 'PCT_TUJ1', 'PCT_S100b', 'PCT_Nestin']
clinical_num = ['ALSFRS_R_Baseline_Value', 'ALSFRS_R_Latest_Value', 'ALSFRS_R_PROGRESSION_SLOPE']
genetic_num = ['C9orf72_repeat_length', 'ATXN2_repeat_length']
all_num_cols = age_cols + ipsc_qc + clinical_num + genetic_num

# Ensure numeric
for col in all_num_cols:
    if col in df_sub.columns:
        df_sub[col] = pd.to_numeric(df_sub[col], errors='coerce')

# Distribution Plot
plt.figure(figsize=(22, 18))
sns.set_style("whitegrid")
for i, col in enumerate(all_num_cols):
    if col not in df_sub.columns: continue
    plt.subplot(4, 4, i + 1)
    sns.histplot(df_sub[col].dropna(), kde=True, color='teal', bins=30)
    plt.title(f'Subject-Level: {col}', fontsize=12)
    plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("$DIST_PLOT")

# Heatmap
plt.figure(figsize=(16, 14))
corr = df_sub[all_num_cols].corr()
sns.heatmap(corr, annot=True, cmap='seismic', center=0, fmt=".2f", linewidths=0.5)
plt.title(f"AnswerALS Subject-Level Correlation (N={len(df_sub)})", fontsize=16)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig("$CORR_PLOT")

# 4. SUBJECT-LEVEL AUDIT TEXT
with open("$COVARIATE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(" AnswerALS Subject-Level Covariate Audit\n")
    f.write(f" UNIQUE PARTICIPANTS (N):      {len(df_sub):<5}\n")
    f.write("============================================================\n\n")

    cat_cols = ['SEX', 'SUBJECT_GROUP', 'SUBJECT_SUBGROUP', 'Site_of_Onset', 'EH_C9orf72', 'has_discordant_sample']
    for col in cat_cols:
        if col not in df_sub.columns: continue
        f.write(f"── {col.upper():<35} {'COUNT':<15} {'PERCENT':<15}\n")
        counts = df_sub[col].value_counts(dropna=False)
        for val, count in counts.items():
            f.write(f"    {str(val):<31} {count:>5} ({ (count/len(df_sub))*100:>4.1f}%)\n")
        f.write("\n")

    f.write("── NUMERICAL COVARIATE STATISTICS (SUBJECT-LEVEL) ──────────\n\n")
    for col in all_num_cols:
        if col in df_sub.columns:
            vals = df_sub[col].dropna()
            f.write(f"{col:<30} (n={len(vals)}, {len(df_sub)-len(vals)} missing)\n")
            f.write(f"    mean: {vals.mean():.2f} | median: {vals.median():.2f} | std: {vals.std():.2f}\n\n")

# 5. TISSUE/OMIC SUMMARY (Still row-level to see total sample availability)
with open("$TISSUE_SUMMARY", "w") as f:
    f.write("AnswerALS Sample Availability (Row-Level Audit)\n\n")
    ct = pd.crosstab(df['Primary_Tissue'].fillna('Unknown'), df['omic'].fillna('None'))
    f.write("TISSUE vs OMIC CROSS-TABULATION:\n")
    f.write(ct.to_string())
EOF

echo "Subject-level summaries and plots generated."
