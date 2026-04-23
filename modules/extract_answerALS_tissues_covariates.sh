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
COVARIATE_SUMMARY="$OUTDIR/subject_level_covariate_summary.txt"
TISSUE_SUMMARY="$OUTDIR/subject_level_tissue_summary.txt"
DIST_PLOT="$OUTDIR/subject_level_distributions.png"
CORR_PLOT="$OUTDIR/subject_level_collinearity_heatmap.png"

# ── STEP 1: SUBJECT-LEVEL AGGREGATION & AUDIT ─────────────────────────────────
echo "Normalizing tissues and generating subject-level multi-omic audit..."

python3 - <<EOF
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# 1. LOAD DATA
df = pd.read_csv("$CONSOLIDATED_META", sep="\t", low_memory=False)

# 2. TISSUE NORMALIZATION
# Standardizing the strings to merge 'NT-cell', 'NT-Cell', 'NT', etc.
def normalize_tissue(t):
    t = str(t).strip()
    if t == 'nan' or not t: return 'Unknown'
    # Merge all variations of Non-T-Cell
    if any(x in t.upper() for x in ['NT-CELL', 'NT_CELL', 'NTCELL', '/NT']):
        return 'PBMC/Non-T-Cell'
    # Merge T-Cell variations
    if 'T-CELL' in t.upper() or 'T_CELL' in t.upper():
        return 'PBMC/T-Cell'
    return t

df['Primary_Tissue'] = df['Primary_Tissue'].apply(normalize_tissue)

# Standardize other key columns
for col in ['omic', 'Participant_ID', 'SEX', 'SUBJECT_GROUP']:
    if col in df.columns:
        df[col] = df[col].astype(str).str.strip().replace('nan', np.nan)

# 3. DEFINE AGGREGATION LOGIC
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
    'ALSFRS_R_Latest_Value': 'min',
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

agg_dict = {k: v for k, v in agg_dict.items() if k in df.columns}
df_sub = df.groupby('Participant_ID').agg(agg_dict).reset_index()

# 4. GENERATE FIGURES (Subject-Level)
age_cols = ['AGE_AT_SYMPTOM_ONSET', 'AGE_AT_DEATH', 'Age_Sample_Collection_Cedars', 'Age_at_First_PBMC_Collection']
ipsc_qc = ['PCT_SMI32', 'PCT_ISL1', 'PCT_NKX61', 'PCT_TUJ1', 'PCT_S100b', 'PCT_Nestin']
clinical_num = ['ALSFRS_R_Baseline_Value', 'ALSFRS_R_Latest_Value', 'ALSFRS_R_PROGRESSION_SLOPE']
genetic_num = ['C9orf72_repeat_length', 'ATXN2_repeat_length']
all_num_cols = age_cols + ipsc_qc + clinical_num + genetic_num

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

# Heatmap (Seismic)
plt.figure(figsize=(16, 14))
corr = df_sub[all_num_cols].corr()
sns.heatmap(corr, annot=True, cmap='seismic', center=0, fmt=".2f", linewidths=0.5)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig("$CORR_PLOT")

# 5. SUBJECT-LEVEL TISSUE/OMIC AUDIT WITH OVERLAP
with open("$TISSUE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(" AnswerALS SUBJECT-LEVEL Omic Availability\n")
    f.write(" (Count of Unique Participants per Tissue)\n")
    f.write("============================================================\n\n")
    
    # Get N per Tissue/Omic
    audit = df.groupby(['Primary_Tissue', 'omic'])['Participant_ID'].nunique().unstack(fill_value=0)
    
    # Calculate G+T Overlap per Tissue using a more robust method
    def get_overlap(group):
        # Find set of unique participants for each omic type in this tissue group
        genomics_subs = set(group[group['omic'] == 'genomics']['Participant_ID'])
        transcriptomics_subs = set(group[group['omic'] == 'transcriptomics']['Participant_ID'])
        # Return the size of the intersection
        return len(genomics_subs.intersection(transcriptomics_subs))

    # Apply the overlap logic across tissues
    overlap_counts = df.groupby('Primary_Tissue').apply(get_overlap, include_groups=False)
    audit['G+T_Overlap'] = overlap_counts
    
    f.write(audit.to_string())
    f.write("\n\n" + "─" * 70 + "\n\n")
    
    # Global Overlap (Subjects with Genomics and Transcriptomics regardless of tissue)
    global_genomics = set(df[df['omic'] == 'genomics']['Participant_ID'])
    global_transcriptomics = set(df[df['omic'] == 'transcriptomics']['Participant_ID'])
    total_overlap = len(global_genomics.intersection(global_transcriptomics))
    
    f.write(f"Global Subjects with BOTH Genomics & Transcriptomics (Any Tissue): {total_overlap}\n")
    f.write(f"Total Unique Participants in dataset: {len(df_sub)}\n")

# 6. SUBJECT-LEVEL COVARIATE SUMMARY
with open("$COVARIATE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(" AnswerALS Subject-Level Covariate Audit\n")
    f.write(f" UNIQUE PARTICIPANTS (N):      {len(df_sub):<5}\n")
    f.write("============================================================\n\n")
    
    cat_cols = ['SEX', 'SUBJECT_GROUP', 'SUBJECT_SUBGROUP', 'Site_of_Onset', 'EH_C9orf72']
    for col in cat_cols:
        if col not in df_sub.columns: continue
        f.write(f"── {col.upper():<35} {'COUNT':<15} {'PERCENT':<15}\n")
        counts = df_sub[col].value_counts(dropna=False)
        for val, count in counts.items():
            f.write(f"    {str(val):<31} {count:>5} ({ (count/len(df_sub))*100:>4.1f}%)\n")
        f.write("\n")
EOF

echo "Process complete. Subject-level audit and figures saved to $OUTDIR"
