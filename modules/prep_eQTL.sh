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
echo "============================================"

# 0. PRE-CLEAN: Remove existing files so awk doesn't report old data if script fails
rm -f $OUTDIR/expression_${TISSUE_DIR}.txt $OUTDIR/snp_${TISSUE_DIR}.txt $OUTDIR/covariates_${TISSUE_DIR}_encoded.txt

##############################################
# STEP 1: Process All Data in Python
##############################################
echo "[Step 1] Aligning Genotypes, Expression, and Covariates..."

python3 << EOF
import pandas as pd
import sys

# 1. LOAD DATA
meta = pd.read_csv("$META", header=None, names=['hra', 'subject', 'tissue'])
meta_tissue = meta[meta['tissue'] == "$TISSUE"].copy()

raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
expr.columns = [c.replace('_', '-') for c in expr.columns]

# 2. FIND OVERLAPPING SUBJECTS (The 275)
# This ignores the 10 samples that have RNA but failed Genotype QC
common_subjects = sorted(list(set(meta_tissue['subject']) & set(raw_df['IID'])))
print(f"DEBUG: Found {len(common_subjects)} common subjects.")

# 3. ALIGN SNPS
snp_final = raw_df[raw_df['IID'].isin(common_subjects)].set_index('IID').drop(columns=['FID']).T
snp_final = snp_final[common_subjects] # Force Column Order

# 4. ALIGN EXPRESSION (Excluding the 10)
# Map HRA -> Subject so we can filter by the 275 Subject IDs
hra_to_sub = dict(zip(meta_tissue['hra'], meta_tissue['subject']))
expr_renamed = expr.rename(columns=hra_to_sub)
# Keep ONLY columns in our 275 common_subjects list
expr_final = expr_renamed[[col for col in common_subjects if col in expr_renamed.columns]]

# 5. ALIGN COVARIATES
cov = pd.read_csv("$COV", sep='\t')
pca = pd.read_csv("$PCA", sep=r'\s+').rename(columns={'#IID': 'IID'})

cov_merged = meta_tissue[meta_tissue['subject'].isin(common_subjects)].merge(cov, left_on='hra', right_on='externalsampleid')
cov_merged = cov_merged.merge(pca, left_on='subject', right_on='IID')

cov_merged['is_als'] = cov_merged['subject_group'].apply(lambda x: 1 if 'ALS' in str(x) else 0)
cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})

cov_cols = ['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
cov_final = cov_merged.set_index('subject')[cov_cols].T
cov_final = cov_final[common_subjects] # Force Column Order

# 6. FINAL VALIDATION & SAVE
print(f"Validation Counts: Expr={expr_final.shape[1]}, SNP={snp_final.shape[1]}, Cov={cov_final.shape[1]}")

if expr_final.shape[1] == snp_final.shape[1] == cov_final.shape[1]:
    expr_final.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')
    snp_final.to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t')
    cov_final.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', quoting=0)
    print(f"SUCCESS: Saved {expr_final.shape[1]} samples.")
else:
    print("FATAL ERROR: Column counts do not match!")
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
