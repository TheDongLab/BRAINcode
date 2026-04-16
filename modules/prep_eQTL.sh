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
# STEP 1: Process All Data in Python
##############################################
echo "[Step 1] Loading, Renaming, and Intersecting all data types..."

python3 << EOF
import pandas as pd
import numpy as np
import sys

# 1. LOAD THE BRIDGE
# meta: RNA ID (hra) -> Subject
# samp_map: DNA ID (hda) -> Subject
meta = pd.read_csv("$META", header=None, names=['hra', 'subject', 'tissue'])
meta = meta[meta['tissue'] == "$TISSUE"]
samp_map = pd.read_csv("$SAMPLE_MAP", sep='\t', header=None, names=['hda', 'subject'])
bridge = pd.merge(meta, samp_map, on='subject')

# 2. PROCESS GENOTYPES
# .raw file uses Short IDs (Subject) in the IID column
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
sub_to_hda = dict(zip(bridge['subject'], bridge['hda']))

# Filter and rename to Long ID (hda)
snp_data = raw_df[raw_df['IID'].isin(bridge['subject'])].copy()
snp_data['IID'] = snp_data['IID'].map(sub_to_hda)
snp_final = snp_data.set_index('IID').drop(columns=['FID']).T

# 3. PROCESS EXPRESSION
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
expr.columns = [c.replace('_', '-') for c in expr.columns]
rna_to_hda = dict(zip(bridge['hra'], bridge['hda']))

# Filter and rename to Long ID (hda)
valid_expr_cols = [c for c in expr.columns if c in rna_to_hda]
expr_final = expr[valid_expr_cols].rename(columns=rna_to_hda)

# 4. PROCESS COVARIATES & PCA
cov = pd.read_csv("$COV", sep='\t')
pca = pd.read_csv("$PCA", sep=r'\s+').rename(columns={'#IID': 'IID'})

# Merge all metadata/PCA onto the bridge
cov_merged = bridge.merge(cov, left_on='hra', right_on='externalsampleid')
pca_to_hda = dict(zip(bridge['subject'], bridge['hda']))
pca['IID_hda'] = pca['IID'].map(pca_to_hda)
cov_merged = cov_merged.merge(pca, left_on='hda', right_on='IID_hda')

# Encodings
def get_als(g): return 1 if 'ALS' in str(g) else 0
cov_merged['is_als'] = cov_merged['subject_group'].apply(get_als)
cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})

cov_cols = ['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
cov_final = cov_merged.set_index('hda')[cov_cols].T

# 5. THE STRICT THREE-WAY INTERSECTION (The Fix)
# Find samples present in all three processed objects
common_samples = sorted(list(set(expr_final.columns) & set(snp_final.columns) & set(cov_final.columns)))

print(f"--- Alignment Results ---")
print(f"Common samples found: {len(common_samples)}")

if len(common_samples) > 0:
    # Save with EXACT same columns in EXACT same order
    expr_final[common_samples].to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')
    snp_final[common_samples].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t')
    cov_final[common_samples].to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', quoting=0)
else:
    print("Error: Intersection failed. Printing ID examples for debug:")
    if not expr_final.empty: print("Expr ID example:", expr_final.columns[0])
    if not snp_final.empty: print("SNP ID example:", snp_final.columns[0])
    if not cov_final.empty: print("Cov ID example:", cov_final.columns[0])
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

##############################################
# STEP 3: Cleanup
##############################################
echo "[Step 3] Cleaning up..."
rm -f $OUTDIR/tmp_*

echo "============================================"
echo "  Prep Complete for $TISSUE"
echo "  $(date)"
echo "============================================"
