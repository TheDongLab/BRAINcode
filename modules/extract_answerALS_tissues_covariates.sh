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

# ── STEP 1: SUBJECT-LEVEL AGGREGATION & AUDIT ─────────────────────────────────
echo "Normalizing DNA Source and calculating Disease Duration..."

python3 - <<EOF
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# 1. LOAD DATA
df = pd.read_csv("$CONSOLIDATED_META", sep="\t", low_memory=False)

# 2. FEATURE ENGINEERING: DISEASE DURATION
# Calculate (Death - Onset) * 12 to get duration in months
for col in ['AGE_AT_DEATH', 'AGE_AT_SYMPTOM_ONSET']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

df['DISEASE_DURATION_IN_MONTHS'] = (df['AGE_AT_DEATH'] - df['AGE_AT_SYMPTOM_ONSET']) * 12

# 3. TISSUE NORMALIZATION (DNA Source)
def normalize_dna_source(s):
    s = str(s).strip()
    if s.lower() in ['nan', 'none', '']: return 'Unknown'
    if any(x in s.upper() for x in ['NT-CELL', 'NT_CELL', 'NTCELL', '/NT']):
        return 'PBMC/Non-T-Cell'
    if 'T-CELL' in s.upper() or 'T_CELL' in s.upper():
        return 'PBMC/T-Cell'
    return s

SOURCE_COL = 'Transcriptomic_Sequencing_DNA_Source'
df['Analyzed_Tissue'] = df[SOURCE_COL].apply(normalize_dna_source) if SOURCE_COL in df.columns else df['Primary_Tissue'].apply(normalize_dna_source)

# 4. DEFINE AGGREGATION LOGIC
# Adding the new line IDs and full repeat lengths
agg_dict = {
    'SEX': 'first',
    'RACE': 'first',
    'SUBJECT_GROUP': 'first',
    'SUBJECT_SUBGROUP': 'first',
    'Site_of_Onset': 'first',
    'iPSCLine': 'first',
    'diMNs_line': 'first',
    'EH_C9orf72': 'max',
    'has_discordant_sample': 'max',
    'AGE_AT_SYMPTOM_ONSET': 'max',
    'AGE_AT_DEATH': 'max',
    'DISEASE_DURATION_IN_MONTHS': 'max',
    'C9orf72_repeat_length': 'max',
    'ATXN2_repeat_length': 'max',
    'C9orf72_repeat_length_full': 'first',
    'ATXN2_repeat_length_full': 'first',
    'ALSFRS_R_PROGRESSION_SLOPE': 'mean',
    'PCT_TUJ1': 'mean',
    'PCT_SMI32': 'mean'
}

agg_dict = {k: v for k, v in agg_dict.items() if k in df.columns}
df_sub = df.groupby('Participant_ID').agg(agg_dict).reset_index()

# 5. GENERATE FIGURES (Subject-Level)
# Include the new duration in plots
num_cols_to_plot = ['AGE_AT_SYMPTOM_ONSET', 'AGE_AT_DEATH', 'DISEASE_DURATION_IN_MONTHS', 
                    'C9orf72_repeat_length', 'ATXN2_repeat_length', 'PCT_TUJ1', 'PCT_SMI32']
num_cols_to_plot = [c for c in num_cols_to_plot if c in df_sub.columns]

for col in num_cols_to_plot:
    df_sub[col] = pd.to_numeric(df_sub[col], errors='coerce')

# Dist Plots
plt.figure(figsize=(20, 15))
sns.set_style("whitegrid")
for i, col in enumerate(num_cols_to_plot):
    plt.subplot(3, 3, i + 1)
    sns.histplot(df_sub[col].dropna(), kde=True, color='teal', bins=30)
    plt.title(f'Subject-Level: {col}', fontsize=12)
plt.tight_layout()
plt.savefig("$DIST_PLOT")

# Heatmap
plt.figure(figsize=(14, 12))
corr = df_sub[num_cols_to_plot].corr()
sns.heatmap(corr, annot=True, cmap='seismic', center=0, fmt=".2f")
plt.title("Subject-Level Clinical/Genetic Correlation")
plt.tight_layout()
plt.savefig("$CORR_PLOT")

# 6. TISSUE SUMMARY (DNA Source with G+T Overlap)
with open("$TISSUE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(" AnswerALS SUBJECT-LEVEL Omic Availability (By DNA Source)\n")
    f.write("============================================================\n\n")
    
    audit = df.groupby(['Analyzed_Tissue', 'omic'])['Participant_ID'].nunique().unstack(fill_value=0)
    
    def get_overlap(group):
        gen_subs = set(group[group['omic'] == 'genomics']['Participant_ID'])
        tra_subs = set(group[group['omic'] == 'transcriptomics']['Participant_ID'])
        return len(gen_subs.intersection(tra_subs))

    audit['G+T_Overlap'] = df.groupby('Analyzed_Tissue').apply(get_overlap, include_groups=False)
    f.write(audit.to_string())

# 7. COVARIATE SUMMARY (With New Variables)
with open("$COVARIATE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(" AnswerALS Subject-Level Covariate Audit\n")
    f.write(f" UNIQUE PARTICIPANTS (N):      {len(df_sub):<5}\n")
    f.write("============================================================\n\n")
    
    cat_cols = ['SEX', 'RACE', 'SUBJECT_GROUP', 'Site_of_Onset', 'iPSCLine', 'diMNs_line']
    for col in cat_cols:
        if col in df_sub.columns:
            f.write(f"── {col.upper():<35}\n")
            counts = df_sub[col].value_counts(dropna=False).head(20) # Head 20 for lines
            for val, count in counts.items():
                f.write(f"    {str(val):<31} {count:>5} ({ (count/len(df_sub))*100:>4.1f}%)\n")
            f.write("\n")

    f.write("── NUMERICAL COVARIATE STATISTICS (SUBJECT-LEVEL) ──────────\n\n")
    for col in num_cols_to_plot:
        vals = df_sub[col].dropna()
        f.write(f"{col:<30} (n={len(vals)}, {len(df_sub)-len(vals)} missing)\n")
        f.write(f"    mean: {vals.mean():.2f} | median: {vals.median():.2f}\n\n")
EOF

echo "Done. Audit and plots generated in $OUTDIR"
