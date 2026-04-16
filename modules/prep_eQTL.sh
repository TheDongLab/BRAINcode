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

# Pre-clean stale files
rm -f $OUTDIR/expression_${TISSUE_DIR}.txt $OUTDIR/snp_${TISSUE_DIR}.txt $OUTDIR/covariates_${TISSUE_DIR}_encoded.txt

python3 << EOF
import pandas as pd
import sys

# 1. LOAD DATA
meta = pd.read_csv("$META", header=None, names=['hra', 'subject', 'tissue'])
meta_tissue = meta[meta['tissue'] == "$TISSUE"].copy()

raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
expr.columns = [c.replace('_', '-') for c in expr.columns]

# 2. DEFINE COMMON SUBJECTS (Dynamic N)
common_subjects = sorted(list(set(meta_tissue['subject']) & set(raw_df['IID'])))
num_common = len(common_subjects)
print(f"DEBUG: Found {num_common} common subjects between Genotypes and Metadata.")

# 3. ALIGN SNPS
snp_final = raw_df[raw_df['IID'].isin(common_subjects)].set_index('IID').drop(columns=['FID']).T
snp_final = snp_final[common_subjects]

# 4. ALIGN EXPRESSION
hra_to_sub = dict(zip(meta_tissue['hra'], meta_tissue['subject']))
valid_hras = [hra for hra, sub in hra_to_sub.items() if sub in common_subjects and hra in expr.columns]
expr_final = expr[valid_hras].rename(columns=hra_to_sub)
expr_final = expr_final[common_subjects]

# 5. ALIGN COVARIATES
cov = pd.read_csv("$COV", sep='\t')
pca = pd.read_csv("$PCA", sep=r'\s+').rename(columns={'#IID': 'IID'})

cov_merged = meta_tissue.merge(cov, left_on='hra', right_on='externalsampleid').merge(pca, left_on='subject', right_on='IID')

cov_merged['is_als'] = cov_merged['subject_group'].apply(lambda x: 1 if 'ALS' in str(x) else 0)
cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})

cov_final = cov_merged[cov_merged['subject'].isin(common_subjects)].set_index('subject')[['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']].T
cov_final = cov_final[common_subjects]

# 6. FINAL VALIDATION (Dynamic Check)
print(f"Validation Counts: Expr={expr_final.shape[1]}, SNP={snp_final.shape[1]}, Cov={cov_final.shape[1]}")

if expr_final.shape[1] == snp_final.shape[1] == cov_final.shape[1] == num_common:
    expr_final.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')
    snp_final.to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t')
    cov_final.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', quoting=0)
    print(f"SUCCESS: Aligned and saved {num_common} samples for $TISSUE.")
else:
    print(f"FATAL ERROR: Column counts do not match the expected count of {num_common}!")
    sys.exit(1)
EOF

### Location Files: ###
echo "Generating location files..."

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
