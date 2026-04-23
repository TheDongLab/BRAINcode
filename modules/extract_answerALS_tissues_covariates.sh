#!/bin/bash
#SBATCH --job-name=extract_answerALS_covariates
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/extract_answerALS_covariates.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/extract_answerALS_covariates.err

set -euo pipefail

# ── PATHS ─────────────────────────────────────────────────────────────────────
DATA_DIR="/home/zw529/donglab/data/answer_ALS"
OUTDIR="$DATA_DIR/QTL"
CONSOLIDATED_META="$DATA_DIR/consolidated_metadata_full.tsv"

mkdir -p "$OUTDIR"

# Output Files
COVARIATE_SUMMARY="$OUTDIR/covariate_summary.txt"
TISSUE_SUMMARY="$OUTDIR/tissue_summary.txt"
DIST_PLOT="$OUTDIR/answerALS_distributions.png"
CORR_PLOT="$OUTDIR/collinearity_heatmap.png"

# ── STEP 1: PROCESSING ────────────────────────────────────────────────────────
echo "Generating split DNA/RNA tissue audits..."

python3 - <<EOF
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# 1. LOAD DATA
df = pd.read_csv("$CONSOLIDATED_META", sep="\t", low_memory=False)

# 2. FEATURE ENGINEERING & NORMALIZATION
for col in ['AGE_AT_DEATH', 'AGE_AT_SYMPTOM_ONSET']:
    df[col] = pd.to_numeric(df[col], errors='coerce')
df['DISEASE_DURATION_IN_MONTHS'] = (df['AGE_AT_DEATH'] - df['AGE_AT_SYMPTOM_ONSET']) * 12

def normalize_tissue(s):
    s = str(s).strip()
    if s.lower() in ['nan', 'none', '']: return 'Unknown'
    if any(x in s.upper() for x in ['NT-CELL', 'NT_CELL', 'NTCELL', '/NT']): return 'PBMC/Non-T-Cell'
    if 'T-CELL' in s.upper() or 'T_CELL' in s.upper(): return 'PBMC/T-Cell'
    return s

# DNA Source (Usually WGS) vs RNA Source (Usually iPSC-diMNs)
df['DNA_Source_Clean'] = df['Primary_Tissue'].apply(normalize_tissue)
df['RNA_Source_Clean'] = df['Transcriptomic_Sequencing_DNA_Source'].apply(normalize_tissue)

# 3. SUBJECT-LEVEL AGGREGATION
agg_dict = {
    'SEX': 'first', 'RACE': 'first', 'SUBJECT_GROUP': 'first', 'Site_of_Onset': 'first',
    'EH_C9orf72': 'max', 'AGE_AT_SYMPTOM_ONSET': 'max', 'AGE_AT_DEATH': 'max',
    'DISEASE_DURATION_IN_MONTHS': 'max', 'Age_Sample_Collection_Cedars': 'max',
    'Age_at_First_PBMC_Collection': 'max', 'ALSFRS_R_Baseline_Value': 'max',
    'ALSFRS_R_Latest_Value': 'min', 'ALSFRS_R_PROGRESSION_SLOPE': 'mean',
    'C9orf72_repeat_length': 'max', 'ATXN2_repeat_length': 'max',
    'PCT_SMI32': 'mean', 'PCT_ISL1': 'mean', 'PCT_NKX61': 'mean',
    'PCT_TUJ1': 'mean', 'PCT_S100b': 'mean', 'PCT_Nestin': 'mean'
}
agg_dict = {k: v for k, v in agg_dict.items() if k in df.columns}
df_sub = df.groupby('Participant_ID').agg(agg_dict).reset_index()

# 4. VARIABLE LIST FOR PLOTS
num_vars = [
    'AGE_AT_SYMPTOM_ONSET', 'AGE_AT_DEATH', 'DISEASE_DURATION_IN_MONTHS',
    'Age_Sample_Collection_Cedars', 'Age_at_First_PBMC_Collection',
    'ALSFRS_R_Baseline_Value', 'ALSFRS_R_Latest_Value', 'ALSFRS_R_PROGRESSION_SLOPE',
    'C9orf72_repeat_length', 'ATXN2_repeat_length',
    'PCT_SMI32', 'PCT_ISL1', 'PCT_NKX61', 'PCT_TUJ1', 'PCT_S100b', 'PCT_Nestin'
]
active_cols = [c for c in num_vars if c in df_sub.columns]
for col in active_cols:
    df_sub[col] = pd.to_numeric(df_sub[col], errors='coerce')

# Distribution Plots
plt.figure(figsize=(26, 22))
sns.set_style("whitegrid")
for i, col in enumerate(active_cols):
    plt.subplot(6, 3, i + 1)
    sns.histplot(df_sub[col].dropna(), kde=True, color='teal', bins=30)
    plt.title(f'Subject-Level: {col}', fontsize=12)
plt.tight_layout()
plt.savefig("$DIST_PLOT")

# Heatmap
plt.figure(figsize=(18, 16))
corr = df_sub[active_cols].corr()
sns.heatmap(corr, annot=True, cmap='seismic', center=0, fmt=".2f", linewidths=0.5)
plt.xticks(rotation=45, ha='right')
plt.savefig("$CORR_PLOT")

# 5. SPLIT TISSUE SUMMARY TABLES
with open("$TISSUE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(" AnswerALS SUBJECT-LEVEL MULTI-OMIC AUDIT\n")
    f.write("============================================================\n\n")

    # Table 1: Genomic Perspective (DNA Source)
    f.write("--- TABLE 1: GENOMIC PERSPECTIVE (Source of DNA/WGS) ---\n")
    dna_audit = df.groupby(['DNA_Source_Clean', 'omic'])['Participant_ID'].nunique().unstack(fill_value=0)
    f.write(dna_audit.to_string())
    f.write("\n\n")

    # Table 2: Transcriptomic Perspective (Source of RNA/iPSC-Lines)
    f.write("--- TABLE 2: TRANSCRIPTOMIC PERSPECTIVE (Source of RNA/Library) ---\n")
    rna_audit = df.groupby(['RNA_Source_Clean', 'omic'])['Participant_ID'].nunique().unstack(fill_value=0)
    
    # Calculate G+T Overlap per RNA tissue (The 'QTL-Ready' count)
    def get_overlap(group):
        gen_subs = set(df[df['omic'] == 'genomics']['Participant_ID'])
        tra_subs = set(group[group['omic'] == 'transcriptomics']['Participant_ID'])
        return len(gen_subs.intersection(tra_subs))

    rna_audit['G+T_Overlap'] = df.groupby('RNA_Source_Clean').apply(get_overlap, include_groups=False)
    f.write(rna_audit.to_string())
    f.write("\n\n")
    
    f.write(f"Global Subject Overlap (G+T Any Tissue): {len(set(df[df['omic']=='genomics']['Participant_ID']).intersection(set(df[df['omic']=='transcriptomics']['Participant_ID'])))}\n")

# 6. COVARIATE SUMMARY
with open("$COVARIATE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(" AnswerALS Subject-Level Covariate Audit\n")
    f.write(f" UNIQUE PARTICIPANTS (N):      {len(df_sub):<5}\n")
    f.write("============================================================\n\n")
    cat_cols = ['SEX', 'RACE', 'SUBJECT_GROUP', 'Site_of_Onset']
    for col in cat_cols:
        if col in df_sub.columns:
            f.write(f"── {col.upper():<35}\n")
            counts = df_sub[col].value_counts(dropna=False)
            for val, count in counts.items():
                f.write(f"    {str(val):<31} {count:>5} ({ (count/len(df_sub))*100:>4.1f}%)\n")
            f.write("\n")
    f.write("── FULL NUMERICAL STATISTICS (SUBJECT-LEVEL) ──────────────\n\n")
    for col in active_cols:
        vals = df_sub[col].dropna()
        f.write(f"{col:<30} (n={len(vals)}, {len(df_sub)-len(vals)} missing)\n")
        f.write(f"    mean: {vals.mean():.2f} | median: {vals.median():.2f} | std: {vals.std():.2f}\n\n")
EOF

echo "Process complete. Split tables and full plots generated."
