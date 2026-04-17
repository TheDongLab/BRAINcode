#!/bin/bash
#SBATCH --job-name=prep_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=399G

set -euo pipefail

TISSUE="$1" # e.g., Motor_Cortex
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')
OUTDIR=/home/zw529/donglab/data/target_ALS/$TISSUE_DIR/eQTL
mkdir -p $OUTDIR

### Paths ###
QTL_DIR=/home/zw529/donglab/data/target_ALS/QTL
PLINK=$QTL_DIR/plink
EXPR=$QTL_DIR/expression_matrix.txt
# Using the two files required for the count logic
METADATA_CSV=/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv
BREAKDOWN_TSV=$QTL_DIR/patient_tissue_breakdown.tsv
COV=$QTL_DIR/covariates.tsv
RAW=$PLINK/joint_autosomes_matrixEQTL.raw
BIM=$PLINK/joint_autosomes_filtered_bed.bim
PCA=$PLINK/joint_pca.eigenvec
GTF_BED6=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6

echo "============================================"
echo "  eQTL Prep for $TISSUE"
echo "============================================"

# Pre-clean stale files
rm -f $OUTDIR/expression_${TISSUE_DIR}.txt $OUTDIR/snp_${TISSUE_DIR}.txt $OUTDIR/covariates_${TISSUE_DIR}_encoded.txt

python3 << EOF
import pandas as pd
import sys
import re

# 1. DEFINE PRECISION PATTERN (TISSUE_REMAP alignment)
patterns = {
    "Motor_Cortex": "Motor Cortex Lateral|Motor Cortex Medial|Lateral Motor Cortex|Medial Motor Cortex|Primary Motor Cortex L|Primary Motor Cortex M|Cortex_Motor_Unspecified|Cortex_Motor_BA4|BA4 Motor Cortex|Lateral_motor_cortex|Motor Cortex|BA4",
    "Cervical_Spinal_Cord": "Spinal_Cord_Cervical|Cervical Spinal Cord|Cervical_spinal_cord|Spinal_cord_Cervical|Cervical",
    "Lumbar_Spinal_Cord": "Lumbar Spinal Cord|Spinal_Cord_Lumbosacral|Lumbosacral_Spinal_Cord|Lumbar_spinal_cord|Lumbar|Lumbosacral",
    "Thoracic_Spinal_Cord": "Thoracic Spinal Cord|Thoracic",
    "Frontal_Cortex": "Frontal Cortex|Frontal",
    "Cerebellum": "Cerebellum"
}
pattern = patterns.get("$TISSUE", "$TISSUE")

# 2. LOAD DATA
# Load Ground Truth Subjects
breakdown = pd.read_csv("$BREAKDOWN_TSV", sep='\t')
# Filter subjects for target tissue
approved_subjects = breakdown[breakdown['tissues'].str.contains("$TISSUE", na=False)]['subject_id'].unique()

# Load and Filter Metadata (Handling column shifts/case issues)
meta = pd.read_csv("$METADATA_CSV", low_memory=False)
# Identify rows: Must be an approved subject AND match the tissue pattern
meta_filt = meta[meta.iloc[:, 1].isin(approved_subjects)].copy() # Column 2 is Subject_ID
# Case-insensitive search across the whole row for the tissue pattern
mask = meta_filt.apply(lambda row: row.astype(str).str.contains(pattern, case=False).any(), axis=1)
meta_tissue = meta_filt[mask].copy()

# Extract RIN (Cols 17/18) and keep best replicate per subject
def get_rin(row):
    try:
        val1 = float(row.iloc[16]) if pd.notnull(row.iloc[16]) else 0
        val2 = float(row.iloc[17]) if pd.notnull(row.iloc[17]) else 0
        return max(val1, val2)
    except: return 0

meta_tissue['RIN_score'] = meta_tissue.apply(get_rin, axis=1)
meta_unique = meta_tissue.sort_values('RIN_score', ascending=False).drop_duplicates(subset=[meta.columns[1]])

# Load Genotypes and Expression
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)

# 3. IDENTIFY TRIPLE INTERSECTION (Subject ID)
# Convert Metadata HRA IDs to match Matrix (CGND-HRA -> CGND_HRA)
meta_unique['hra_clean'] = meta_unique.iloc[:, 0].str.replace('-', '_')

# We need Subject IDs present in breakdown, metadata, and genotype IIDs
common_subjects = sorted(list(set(meta_unique.iloc[:, 1]) & set(raw_df['IID'])))

# Final alignment check: Metadata HRA must exist in Matrix columns (including _x/_y)
sub_to_hra = dict(zip(meta_unique.iloc[:, 1], meta_unique['hra_clean']))
final_aligned_subjects = []
for s in common_subjects:
    hra = sub_to_hra[s]
    # Check for direct match or technical replicates in expr matrix
    matches = [c for c in expr.columns if c.startswith(hra)]
    if matches:
        final_aligned_subjects.append(s)
        sub_to_hra[s] = matches[0] # Map subject to specific matrix column

num_final = len(final_aligned_subjects)

# 4. ALIGN FILES
# Expression
expr_final = expr[ [sub_to_hra[s] for s in final_aligned_subjects] ]
expr_final.columns = final_aligned_subjects

# SNPs
snp_final = raw_df[raw_df['IID'].isin(final_aligned_subjects)].set_index('IID').drop(columns=['FID']).T
def clean_snp_id(full_id):
    parts = str(full_id).split(':')
    return f"chr{parts[0]}:{parts[1]}" if len(parts) >= 2 else full_id
snp_final.index = [clean_snp_id(idx) for idx in snp_final.index]

# Covariates
cov = pd.read_csv("$COV", sep='\t')
pca = pd.read_csv("$PCA", sep=r'\s+').rename(columns={'#IID': 'IID'})

# Merge metadata with covariates (on HRA) and PCA (on Subject)
meta_aligned = meta_unique[meta_unique.iloc[:, 1].isin(final_aligned_subjects)]
cov_merged = meta_aligned.merge(cov, left_on=meta.columns[0], right_on='externalsampleid').merge(pca, left_on=meta.columns[1], right_on='IID')

cov_merged['is_als'] = cov_merged['subject_group'].apply(lambda x: 1 if 'ALS' in str(x) else 0)
cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})

cov_final = cov_merged.set_index(meta.columns[1]).reindex(final_aligned_subjects)
cov_final = cov_final[['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']].T

# 5. SAVE
print(f"Validation Counts: Expr={expr_final.shape[1]}, SNP={snp_final[final_aligned_subjects].shape[1]}, Cov={cov_final.shape[1]}")

if expr_final.shape[1] == num_final and snp_final.shape[1] == num_final:
    expr_final.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="geneid")
    snp_final[final_aligned_subjects].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="snpid")
    cov_final.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', index=True, index_label="id", quoting=0)
    print(f"SUCCESS: Saved {num_final} aligned samples for $TISSUE.")
else:
    print("FATAL ERROR: Mismatch detected!")
    sys.exit(1)
EOF

### Location Files ###
echo "Generating location files..."
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt

echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3}' $GTF_BED6 >> $OUTDIR/gene_location.txt

echo "============================================"
echo "  Prep Complete for $TISSUE at $(date)"
echo "============================================"
