#!/bin/bash
#SBATCH --job-name=prep_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=399G

set -euo pipefail

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')
OUTDIR=/home/zw529/donglab/data/target_ALS/$TISSUE_DIR/eQTL
mkdir -p $OUTDIR

### Paths ###
QTL_DIR=/home/zw529/donglab/data/target_ALS/QTL
PLINK=$QTL_DIR/plink
EXPR=$QTL_DIR/expression_matrix.txt
META=$QTL_DIR/expression_sample_metadata.csv
COV=$QTL_DIR/covariates.tsv
RAW=$PLINK/joint_autosomes_matrixEQTL.raw
BIM=$PLINK/joint_autosomes_filtered_bed.bim
PCA=$PLINK/joint_pca.eigenvec
GTF_BED6=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6

echo "============================================"
echo "  eQTL Prep for $TISSUE"
echo "  Targeting Subject-Level Alignment"
echo "============================================"

##############################################
# STEP 1: Process All Data in Python
##############################################
echo "[Step 1] Aligning Genotypes, Expression, and Covariates..."

python3 << EOF
import pandas as pd
import sys

# 1. LOAD METADATA
# This is our master key: Column 0 is the RNA ID, Column 1 is the Subject ID
meta = pd.read_csv("$META", header=None, names=['hra', 'subject', 'tissue'])
meta_tissue = meta[meta['tissue'] == "$TISSUE"].copy()

# 2. LOAD GENOTYPES (The Anchor)
# The .raw file uses Subject IDs (NEU..., JHU..., SD...) in the 'IID' column
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])

# 3. IDENTIFY COMMON SUBJECTS
# Find subjects who have BOTH genotype data AND tissue-specific RNA-seq data
common_subjects = set(meta_tissue['subject']) & set(raw_df['IID'])
print(f"DEBUG: Found {len(common_subjects)} subjects with both Genotypes and $TISSUE Expression.")

# 4. PROCESS SNPS
# Filter genotypes to common subjects and set IID (Subject ID) as index
snp_final = raw_df[raw_df['IID'].isin(common_subjects)].copy()
snp_final = snp_final.set_index('IID').drop(columns=['FID']).T

# 5. PROCESS EXPRESSION
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
expr.columns = [c.replace('_', '-') for c in expr.columns]

# Create a mapping from the RNA sample ID to the Subject ID
hra_to_sub = dict(zip(meta_tissue['hra'], meta_tissue['subject']))

# Select columns that are in our common subject list
# We rename them to Subject IDs so they match the SNP file exactly
valid_expr_cols = [c for c in expr.columns if c in hra_to_sub and hra_to_sub[c] in common_subjects]
expr_final = expr[valid_expr_cols].rename(columns=hra_to_sub)

# 6. PROCESS COVARIATES & PCA
cov = pd.read_csv("$COV", sep='\t')
pca = pd.read_csv("$PCA", sep=r'\s+').rename(columns={'#IID': 'IID'})

# Merge covariates and PCA using Subject ID as the join key
cov_merged = meta_tissue[meta_tissue['subject'].isin(common_subjects)].merge(cov, left_on='hra', right_on='externalsampleid')
cov_merged = cov_merged.merge(pca, left_on='subject', right_on='IID')

# Binary encodings
cov_merged['is_als'] = cov_merged['subject_group'].apply(lambda x: 1 if 'ALS' in str(x) else 0)
cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})

cov_cols = ['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
cov_final = cov_merged.set_index('subject')[cov_cols].T

# 7. FINAL ALIGNMENT AND SAVE
# We use the sorted list of common subjects to ensure Column order is identical
final_ids = sorted(list(common_subjects))

print(f"DEBUG: Final Aligned Sample Count: {len(final_ids)}")

if len(final_ids) > 0:
    expr_final[final_ids].to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')
    snp_final[final_ids].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t')
    cov_final[final_ids].to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', quoting=0)
else:
    print("FATAL: No overlapping samples found. Check metadata and raw file IDs.")
    sys.exit(1)
EOF

##############################################
# STEP 2: Location Files
##############################################
echo "[Step 2] Generating location files..."

# SNP Locations
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt

# Gene Locations
echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3}' $GTF_BED6 >> $OUTDIR/gene_location.txt

echo "============================================"
echo "  Prep Complete for $TISSUE"
echo "  $(date)"
echo "============================================"
