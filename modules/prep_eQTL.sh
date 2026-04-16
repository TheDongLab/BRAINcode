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
SAMPLE_MAP=$QTL_DIR/sample_map.txt
COV=$QTL_DIR/covariates.tsv
RAW=$PLINK/joint_autosomes_matrixEQTL.raw
BIM=$PLINK/joint_autosomes_filtered_bed.bim
PCA=$PLINK/joint_pca.eigenvec
GTF_BED6=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6

echo "============================================"
echo "  eQTL Prep for $TISSUE"
echo "  $(date)"
echo "============================================"

##############################################
# STEP 1: Pre-process Covariates
##############################################
echo "[Step 1] Cleaning covariates..."

python3 << EOF
import pandas as pd
df = pd.read_csv("$COV", sep='\t')
# Keep only those with necessary metadata
df = df[df['externalsampleid'].notna()]

def collapse(g):
    if 'ALS' in str(g): return 1
    if 'Non-Neurological' in str(g): return 0
    return None

df['is_als'] = df['subject_group'].apply(collapse)
df = df[df['is_als'].notna()]
df.to_csv("$OUTDIR/tmp_clean_cov.tsv", sep='\t', index=False)
EOF

##############################################
# STEP 2: ID Alignment (DNA <-> RNA)
##############################################
echo "[Step 2] Bridging IDs..."

python3 << EOF
import pandas as pd
import sys

# Load Metadata
meta_df = pd.read_csv("$META", header=None, names=['hra', 'subject', 'tissue'])
meta_filt = meta_df[meta_df['tissue'] == "$TISSUE"].copy()

# Load Sample Map
map_df = pd.read_csv("$SAMPLE_MAP", sep='\t', header=None, names=['hda', 'subject'])

# Bridge them
triplets = pd.merge(meta_filt, map_df, on='subject')
triplets.to_csv("$OUTDIR/tmp_triplets.csv", index=False)
EOF

##############################################
# STEP 3: Initial Data Intersection
##############################################
echo "[Step 3] Initial data loading..."

python3 << EOF
import pandas as pd
import numpy as np

triplets = pd.read_csv("$OUTDIR/tmp_triplets.csv")
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
pca_raw = pd.read_csv("$PCA", sep=r'\s+').rename(columns={'#IID': 'IID'})

# 1. Process SNPs
# Map DNA_ID (hda) to Subject for SNP file
id_map = dict(zip(triplets['hda'], triplets['subject']))
snp_out = raw_df[raw_df['IID'].isin(triplets['hda'])].copy()
snp_out['SUB_ID'] = snp_out['IID'].map(id_map)
snp_out = snp_out.set_index('SUB_ID').drop(columns=['FID', 'IID']).T
snp_out.to_csv("$OUTDIR/tmp_snp_raw.txt", sep='\t')

# 2. Process Expression
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
expr.columns = [c.replace('_', '-') for c in expr.columns]
# Standardize Expr headers to Subject
rna_to_sub = dict(zip(triplets['hra'], triplets['subject']))
expr_subset = expr[[c for c in triplets['hra'] if c in expr.columns]].copy()
expr_subset.rename(columns=rna_to_sub, inplace=True)
expr_subset.to_csv("$OUTDIR/tmp_expr_raw.txt", sep='\t')

# 3. Process Covariates
cov = pd.read_csv("$OUTDIR/tmp_clean_cov.tsv", sep='\t')
cov_final = triplets.merge(cov, left_on='hra', right_on='externalsampleid')
cov_final = cov_final.merge(pca_raw, left_on='subject', right_on='IID')
cov_final['sex_bin'] = cov_final['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})
cov_cols = ['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
cov_out = cov_final.set_index('subject')[cov_cols].T
cov_out.to_csv("$OUTDIR/tmp_cov_raw.txt", sep='\t')
EOF

##############################################
# STEP 4: STRICT FINAL ALIGNMENT
##############################################
echo "[Step 4] Performing strict three-way sample alignment..."

python3 << EOF
import pandas as pd

snp = pd.read_csv("$OUTDIR/tmp_snp_raw.txt", sep='\t', index_col=0)
expr = pd.read_csv("$OUTDIR/tmp_expr_raw.txt", sep='\t', index_col=0)
cov = pd.read_csv("$OUTDIR/tmp_cov_raw.txt", sep='\t', index_col=0)

# Find samples present in all three
common = sorted(list(set(snp.columns) & set(expr.columns) & set(cov.columns)))

print(f"  Final Aligned Sample Count: {len(common)}")

# Save final files with identical column order
snp[common].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t')
expr[common].to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')
cov[common].to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t')
EOF

##############################################
# STEP 5: Location Files
##############################################
echo "[Step 5] Generating location files..."
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt

echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3}' $GTF_BED6 >> $OUTDIR/gene_location.txt

##############################################
# STEP 6: Cleanup
##############################################
echo "[Step 6] Cleaning up..."
rm -f $OUTDIR/tmp_*

echo "============================================"
echo "  Prep Complete for $TISSUE"
echo "  Final count verified."
echo "============================================"
