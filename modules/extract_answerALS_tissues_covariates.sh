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

# ── STEP 1: CONSOLIDATED PROCESSING ──────────────────────────────────────────
echo "Running full subject-level audit and generating distribution plots..."

python3 - <<EOF
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

# 1. LOAD DATA
df = pd.read_csv("$CONSOLIDATED_META", sep="\t", low_memory=False)

# 2. FEATURE ENGINEERING & NORMALIZATION
for col in ['AGE_AT_DEATH', 'AGE_AT_SYMPTOM_ONSET']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Calculate Disease Duration in Months
df['DISEASE_DURATION_IN_MONTHS'] = (df['AGE_AT_DEATH'] - df['AGE_AT_SYMPTOM_ONSET']) * 12

def normalize_tissue(s):
    s = str(s).strip()
    if s.lower() in ['nan', 'none', '']: return 'Unknown'
    if any(x in s.upper() for x in ['NT-CELL', 'NT_CELL', 'NTCELL', '/NT']): return 'PBMC/Non-T-Cell'
    if 'T-CELL' in s.upper() or 'T_CELL' in s.upper(): return 'PBMC/T-Cell'
    return s

df['DNA_Source_Clean'] = df['Primary_Tissue'].apply(normalize_tissue)
df['RNA_Source_Clean'] = df['Transcriptomic_Sequencing_DNA_Source'].apply(normalize_tissue)

# 3. DEFINE NUMERICAL VARIABLES (Fixed the missing definition)
all_num_vars = [
    'AGE_AT_SYMPTOM_ONSET', 'AGE_AT_DEATH', 'DISEASE_DURATION_IN_MONTHS',
    'Age_Sample_Collection_Cedars', 'Age_at_First_PBMC_Collection',
    'ALSFRS_R_Baseline_Value', 'ALSFRS_R_Latest_Value', 'ALSFRS_R_PROGRESSION_SLOPE',
    'C9orf72_repeat_length', 'ATXN2_repeat_length',
    'PCT_SMI32', 'PCT_ISL1', 'PCT_NKX61', 'PCT_TUJ1', 'PCT_S100b', 'PCT_Nestin'
]

# 4. SUBJECT-LEVEL AGGREGATION
agg_dict = {
    'SEX': 'first', 'RACE': 'first', 'SUBJECT_GROUP': 'first', 'Site_of_Onset': 'first',
    'EH_C9orf72': 'max', 'C9orf72_repeat_length_full': 'first', 'ATXN2_repeat_length_full': 'first'
}

# Add numerical vars to aggregation dictionary
for col in all_num_vars:
    if col in df.columns:
        if 'PCT' in col or 'SLOPE' in col:
            agg_dict[col] = 'mean' # Average QC metrics for the subject
        else:
            agg_dict[col] = 'max' # Take max for clinical/age/repeat values

agg_dict = {k: v for k, v in agg_dict.items() if k in df.columns}
df_sub = df.groupby('Participant_ID').agg(agg_dict).reset_index()

# 5. GENERATE DISTRIBUTION FIGURES
active_cols = [c for c in all_num_vars if c in df_sub.columns]
for col in active_cols:
    df_sub[col] = pd.to_numeric(df_sub[col], errors='coerce')

fig = plt.figure(figsize=(26, 26))
sns.set_style("whitegrid")
gs = GridSpec(6, 3, figure=fig)

for i, col in enumerate(active_cols):
    row, loc = divmod(i, 3)
    
    if col == 'C9orf72_repeat_length':
        # Broken Axis subgrid for C9 outliers
        sub_gs = gs[row, loc].subgridspec(1, 2, width_ratios=[3, 1], wspace=0.1)
        ax_main = fig.add_subplot(sub_gs[0])
        ax_outlier = fig.add_subplot(sub_gs[1], sharey=ax_main)
        
        data = df_sub[col].dropna()
        sns.histplot(data, bins=20, binrange=(0, 40), ax=ax_main, color='teal', kde=False)
        sns.histplot(data, bins=20, binrange=(250, 1250), ax=ax_outlier, color='red', alpha=0.6)
        
        ax_main.set_xlim(0, 50)
        ax_outlier.set_xlim(400, 1250)
        ax_main.set_title("C9 Zoom (0-50)")
        ax_outlier.set_title("Outliers")
        
        ax_main.spines['right'].set_visible(False)
        ax_outlier.spines['left'].set_visible(False)
        ax_outlier.yaxis.set_visible(False)
    else:
        ax = fig.add_subplot(gs[row, loc])
        sns.histplot(df_sub[col].dropna(), kde=True, color='teal', bins=30, ax=ax)
        ax.set_title(f'Subject-Level: {col}', fontsize=12)

plt.tight_layout()
plt.savefig("$DIST_PLOT")

# 6. GENERATE CORRELATION HEATMAP
plt.figure(figsize=(18, 16))
corr = df_sub[active_cols].corr()
sns.heatmap(corr, annot=True, cmap='seismic', center=0, fmt=".2f", linewidths=0.5)
plt.title(f"AnswerALS Global Subject-Level Correlation (N={len(df_sub)})", fontsize=16)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig("$CORR_PLOT")

# 7. TISSUE SUMMARY (Split DNA/RNA tables)
with open("$TISSUE_SUMMARY", "w") as f:
    f.write("============================================================\n")
    f.write(" AnswerALS SUBJECT-LEVEL MULTI-OMIC AUDIT\n")
    f.write("============================================================\n\n")

    f.write("--- TABLE 1: GENOMIC PERSPECTIVE (Source of DNA/WGS) ---\n")
    dna_audit = df.groupby(['DNA_Source_Clean', 'omic'])['Participant_ID'].nunique().unstack(fill_value=0)
    f.write(dna_audit.to_string())
    f.write("\n\n")

    f.write("--- TABLE 2: TRANSCRIPTOMIC PERSPECTIVE (Source of RNA/Library) ---\n")
    rna_audit = df.groupby(['RNA_Source_Clean', 'omic'])['Participant_ID'].nunique().unstack(fill_value=0)
    
    def get_overlap(group):
        gen_subs = set(df[df['omic'] == 'genomics']['Participant_ID'])
        tra_subs = set(group[group['omic'] == 'transcriptomics']['Participant_ID'])
        return len(gen_subs.intersection(tra_subs))

    rna_audit['G+T_Overlap'] = df.groupby('RNA_Source_Clean').apply(get_overlap, include_groups=False)
    f.write(rna_audit.to_string())
    f.write("\n\n")
    
    global_gen = set(df[df['omic']=='genomics']['Participant_ID'])
    global_tra = set(df[df['omic']=='transcriptomics']['Participant_ID'])
    f.write(f"Global Subject Overlap (G+T Any Tissue): {len(global_gen.intersection(global_tra))}\n")

# 8. COVARIATE SUMMARY
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

echo "Done. All files generated in $OUTDIR"
