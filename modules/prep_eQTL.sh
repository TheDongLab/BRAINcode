#!/bin/bash
#SBATCH --job-name=prep_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=499G

set -euo pipefail

### Arguments ###
if [ $# -lt 1 ]; then
    echo "ERROR: Usage: sbatch prep_eQTL.sh \"Tissue Name\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

OUTDIR=/home/zw529/donglab/data/target_ALS/$TISSUE_DIR/eQTL
mkdir -p $OUTDIR

### Paths ###
BASE=/home/zw529/donglab/data/target_ALS
QTL_DIR=$BASE/QTL
PLINK=$QTL_DIR/plink
REFS=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode

EXPR=$QTL_DIR/expression_matrix.txt
META=$QTL_DIR/expression_sample_metadata.csv
SAMPLE_MAP=$QTL_DIR/sample_map.txt        
COV=$QTL_DIR/covariates.tsv               
RAW=$PLINK/joint_autosomes_matrixEQTL.raw
BIM=$PLINK/joint_autosomes_filtered_bed.bim
PCA=$PLINK/joint_pca.eigenvec
GTF_BED6=$REFS/gencode.v49.annotation.gene.bed6

echo "============================================"
echo "  eQTL Prep for $TISSUE"
echo "  $(date)"
echo "============================================"

##############################################
# STEP 1: Build ID Mapping (DNA <-> RNA)
##############################################
echo "[1] Filtering metadata..."

python3 << EOF
import pandas as pd
df = pd.read_csv("$COV", sep='\t')
bad_samples = df[df['rin'].isna() & (df['rna_skew'] > 1.0)]['externalsampleid'].tolist()
df = df[~df['externalsampleid'].isin(bad_samples)]

def collapse(g):
    if 'ALS' in str(g): return 1
    if 'Non-Neurological' in str(g): return 0
    return None

df['is_als'] = df['subject_group'].apply(collapse)
df = df[df['is_als'].notna()]
df.to_csv("$OUTDIR/tmp_clean_cov.tsv", sep='\t', index=False)
EOF

sed 's/[[:space:]]\+/,/g' $SAMPLE_MAP | sort -t',' -k2,2 > $OUTDIR/tmp_map.csv
tail -n +2 $META | awk -F',' -v tissue="$TISSUE" '$3 == tissue {print $1","$2}' | sort -t',' -k1,1 > $OUTDIR/tmp_meta.csv
join -t',' -1 1 -2 2 $OUTDIR/tmp_meta.csv $OUTDIR/tmp_map.csv > $OUTDIR/tmp_triplets.csv

##############################################
# STEP 2-6: Alignment and Matrix Generation
##############################################
echo "[2-6] Aligning SNPs, Expression, and Covariates..."

python3 << EOF
import pandas as pd
import numpy as np
import re
import sys

mapping = pd.read_csv("$OUTDIR/tmp_triplets.csv", header=None, names=['sub','hra','hda'])

# FIX 1: Use raw string r'\s+' and cast to string before .str accessor
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
raw_df['FID'] = raw_df['FID'].astype(str).str.replace(r'-b\d+$', '', regex=True)

expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
expr.columns = [c.replace('_', '-') for c in expr.columns]

cov = pd.read_csv("$OUTDIR/tmp_clean_cov.tsv", sep='\t')

# FIX 2: Repeat raw string and string cast for PCA
pca = pd.read_csv("$PCA", sep=r'\s+', header=None).iloc[:, [1,2,3,4,5,6]]
pca.columns = ['ID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
pca['ID'] = pca['ID'].astype(str).str.replace(r'-b\d+$', '', regex=True)

overlap = mapping[
    mapping['hda'].astype(str).isin(raw_df['FID']) & 
    mapping['hra'].isin(expr.columns) &
    mapping['hra'].isin(cov['externalsampleid']) &
    mapping['hda'].astype(str).isin(pca['ID'])
].copy()

print(f"  Final sample overlap count: {len(overlap)}")

if len(overlap) == 0:
    print("ERROR: Zero samples overlapping. Check ID formats.")
    sys.exit(1)

expr_aligned = expr[overlap['hra']]
expr_aligned.columns = overlap['hda']
expr_aligned.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')

snp_aligned = raw_df[raw_df['FID'].isin(overlap['hda'])].set_index('FID').drop(columns=['IID']).T
af = np.nanmean(snp_aligned.values, axis=1) / 2.0
maf = np.minimum(af, 1 - af)
keep_snps = np.where(np.nan_to_num(maf) >= 0.05)[0]
snp_aligned.iloc[keep_snps, :].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t')

cov_final = overlap.merge(cov, left_on='hra', right_on='externalsampleid')
cov_final = cov_final.merge(pca, left_on='hda', right_on='ID')
cov_final['sex_binary'] = cov_final['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})

cov_cols = ['sex_binary', 'age_at_death', 'rin', 'rna_skew', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
cov_final.set_index('hda')[cov_cols].T.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t')
EOF

##############################################
# STEP 7-8: Location Files
##############################################
echo "[7-8] Generating location files..."
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt

echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {
    gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3
}' $GTF_BED6 >> $OUTDIR/gene_location.txt

rm -f $OUTDIR/tmp_*
echo "============================================"
echo "  Prep complete: $OUTDIR"
echo "  $(date)"
echo "============================================"

##############################################
# STEP 9-10: Cleanup and Summary
##############################################
echo "[9] Cleaning up temp files..."
rm -f $OUTDIR/tmp_*

echo "============================================"
echo "  Prep complete for : $TISSUE"
echo "  Output directory  : $OUTDIR"
echo ""
echo "  snp_${TISSUE_DIR}.txt                  -> SNP_file_name"
echo "  expression_${TISSUE_DIR}.txt           -> expression_file_name"
echo "  covariates_${TISSUE_DIR}_encoded.txt   -> covariates_file_name"
echo "  gene_location.txt                      -> gene_location_file_name"
echo "  snp_location.txt                       -> snp_location_file_name"
echo "  $(date)"
echo "============================================"
