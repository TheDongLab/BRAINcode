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

# Writing directly to the tissue's eQTL directory (flattened as requested)
OUTDIR=/home/zw529/donglab/data/target_ALS/$TISSUE_DIR/eQTL
mkdir -p $OUTDIR

### Verified Paths from cluster inspection ###
BASE=/home/zw529/donglab/data/target_ALS
QTL_DIR=$BASE/QTL
PLINK=$QTL_DIR/plink
REFS=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode

# Input Files
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
echo "[1] Filtering metadata and normalizing ID maps..."

# Process Covariates and create binary filters
python3 << EOF
import pandas as pd
df = pd.read_csv("$COV", sep='\t')

# Quality Filter: Drop samples with low quality metrics
bad_samples = df[df['rin'].isna() & (df['rna_skew'] > 1.0)]['externalsampleid'].tolist()
df = df[~df['externalsampleid'].isin(bad_samples)]

# Case/Control Encoding (Matrix eQTL needs numeric)
def collapse(g):
    if 'ALS' in str(g): return 1
    if 'Non-Neurological' in str(g): return 0
    return None

df['is_als'] = df['subject_group'].apply(collapse)
df = df[df['is_als'].notna()]
df.to_csv("$OUTDIR/tmp_clean_cov.tsv", sep='\t', index=False)
EOF

# Normalize IDs across files:
# 1. Map: DNA_ID, Subject_ID (handling the spaces/tabs from your ls output)
sed 's/[[:space:]]\+/,/g' $SAMPLE_MAP | sort -t',' -k2,2 > $OUTDIR/tmp_map.csv

# 2. Metadata: Subject_ID, RNA_ID (HRA)
tail -n +2 $META | awk -F',' -v tissue="$TISSUE" '$3 == tissue {print $1","$2}' | sort -t',' -k1,1 > $OUTDIR/tmp_meta.csv

# 3. Join to create: Subject, RNA_ID, DNA_ID
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

# Load Mapping
mapping = pd.read_csv("$OUTDIR/tmp_triplets.csv", header=None, names=['sub','hra','hda'])

# 1. Load Genotypes
# Using the FID (Col 0) for DNA IDs per your cluster inspection
raw_df = pd.read_csv("$RAW", sep='\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
# Strip DNA versioning (e.g., -b38) to match mapping
raw_df['FID'] = raw_df['FID'].str.replace(r'-b\d+$', '', regex=True)

# 2. Load Expression
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
# Normalize RNA IDs: Matrix uses underscores (CGND_HRA), Mapping uses hyphens (CGND-HRA)
expr.columns = [c.replace('_', '-') for c in expr.columns]

# 3. Load Covariates and PCA
cov = pd.read_csv("$OUTDIR/tmp_clean_cov.tsv", sep='\t')
pca = pd.read_csv("$PCA", sep='\s+', header=None).iloc[:, [1,2,3,4,5,6]]
pca.columns = ['ID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
pca['ID'] = pca['ID'].str.replace(r'-b\d+$', '', regex=True)

# 4. Find the Overlap (Intersection of all 4 data sources)
overlap = mapping[
    mapping['hda'].isin(raw_df['FID']) & 
    mapping['hra'].isin(expr.columns) &
    mapping['hra'].isin(cov['externalsampleid']) &
    mapping['hda'].isin(pca['ID'])
].copy()

print(f"  Final sample overlap count: {len(overlap)}")

if len(overlap) == 0:
    print("ERROR: Zero samples overlapping. Check ID formats.")
    sys.exit(1)

# 5. Export Aligned Expression (Column names = DNA IDs for Matrix eQTL)
expr_aligned = expr[overlap['hra']]
expr_aligned.columns = overlap['hda']
expr_aligned.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')

# 6. Export Aligned SNPs
snp_aligned = raw_df[raw_df['FID'].isin(overlap['hda'])].set_index('FID').drop(columns=['IID']).T
# Filter for MAF > 0.05
af = np.nanmean(snp_aligned.values, axis=1) / 2.0
maf = np.minimum(af, 1 - af)
keep_snps = np.where(np.nan_to_num(maf) >= 0.05)[0]
snp_aligned.iloc[keep_snps, :].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t')

# 7. Export Aligned Covariates
cov_final = overlap.merge(cov, left_on='hra', right_on='externalsampleid')
cov_final = cov_final.merge(pca, left_on='hda', right_on='ID')
cov_final['sex_binary'] = cov_final['sex'].str.lower().map({'male': 1, 'female': 0})

# Transpose for Matrix eQTL
cov_cols = ['sex_binary', 'age_at_death', 'rin', 'rna_skew', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
cov_final.set_index('hda')[cov_cols].T.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t')
EOF

##############################################
# STEP 7-8: Location Files & Validation
##############################################
echo "[7] Generating SNP location file..."
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt

echo "[8] Generating gene location file..."
echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {
    gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3
}' $GTF_BED6 >> $OUTDIR/gene_location.txt

# Final ID overlap validation check
python3 << EOF
expr_file = "$OUTDIR/expression_${TISSUE_DIR}.txt"
loc_file  = "$OUTDIR/gene_location.txt"
with open(loc_file) as f:
    f.readline()
    loc_genes = set(line.split('\t')[0] for line in f)
with open(expr_file) as f:
    f.readline()
    expr_genes = [line.split('\t')[0] for line in f]
found = sum(1 for g in expr_genes if g in loc_genes)
print(f"  Validation: Found location for {found} / {len(expr_genes)} expression genes.")
EOF

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
