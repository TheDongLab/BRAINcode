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

# Critical Step: Ensure 1-to-1 Subject-to-HRA mapping for this tissue
meta_unique = meta_tissue.drop_duplicates(subset=['subject'], keep='first')

raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
expr.columns = [c.replace('_', '-') for c in expr.columns]

# 2. IDENTIFY THE TRUE INTERSECTION
common_subjects = sorted(list(set(meta_unique['subject']) & set(raw_df['IID'])))
num_common = len(common_subjects)
print(f"DEBUG: True intersection count = {num_common}")

# 3. ALIGN SNPS & CLEAN IDs
snp_final = raw_df[raw_df['IID'].isin(common_subjects)].set_index('IID').drop(columns=['FID']).T

# Function to fix IDs: 1:788418:CAG:C_CAG -> chr1:788418
def clean_snp_id(full_id):
    parts = str(full_id).split(':')
    if len(parts) >= 2:
        return f"chr{parts[0]}:{parts[1]}"
    return full_id

snp_final.index = [clean_snp_id(idx) for idx in snp_final.index]

# 4. ALIGN EXPRESSION
sub_to_hra = dict(zip(meta_unique['subject'], meta_unique['hra']))
final_aligned_subjects = [s for s in common_subjects if sub_to_hra[s] in expr.columns]
num_final = len(final_aligned_subjects)

expr_final = expr[[sub_to_hra[s] for s in final_aligned_subjects]]
expr_final.columns = final_aligned_subjects

# 5. ALIGN COVARIATES
cov = pd.read_csv("$COV", sep='\t')
pca = pd.read_csv("$PCA", sep=r'\s+').rename(columns={'#IID': 'IID'})

cov_merged = meta_unique.merge(cov, left_on='hra', right_on='externalsampleid').merge(pca, left_on='subject', right_on='IID')
cov_merged['is_als'] = cov_merged['subject_group'].apply(lambda x: 1 if 'ALS' in str(x) else 0)
cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})

cov_final = cov_merged.set_index('subject').reindex(final_aligned_subjects)
cov_final = cov_final[['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']].T

# 6. FINAL TRIPLE-CHECK & SAVE WITH LABELS
print(f"Validation Counts: Expr={expr_final.shape[1]}, SNP={snp_final[final_aligned_subjects].shape[1]}, Cov={cov_final.shape[1]}")

if expr_final.shape[1] == snp_final[final_aligned_subjects].shape[1] == cov_final.shape[1] == num_final:
    # Adding index=True and index_label ensures Matrix eQTL matches data to locations
    expr_final.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="geneid")
    snp_final[final_aligned_subjects].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="snpid")
    cov_final.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', index=True, index_label="id", quoting=0)
    print(f"SUCCESS: Saved {num_final} aligned samples for $TISSUE.")
else:
    print("FATAL ERROR: Mismatch detected in final reconstruction!")
    sys.exit(1)
EOF

### Location Files ###
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
