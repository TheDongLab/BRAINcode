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

##############################################
# STEP 2: ID Alignment (DNA <-> RNA)
##############################################
echo "[Step 2] Bridging IDs in Python..."

python3 << EOF
import pandas as pd
import sys

# Load Metadata (RNA ID, Subject ID, Tissue)
meta_df = pd.read_csv("$META", header=None, names=['hra', 'subject', 'tissue'])
meta_filt = meta_df[meta_df['tissue'] == "$TISSUE"].copy()

# Load Sample Map (DNA ID, Subject ID)
map_df = pd.read_csv("$SAMPLE_MAP", sep='\t', header=None, names=['hda', 'subject'])

# Bridge them
triplets = pd.merge(meta_filt, map_df, on='subject')
triplets.to_csv("$OUTDIR/tmp_triplets.csv", index=False)
print(f"  Mapped {len(triplets)} candidates for $TISSUE.")
EOF

##############################################
# STEP 3: Load and Intersect High-D Data
##############################################
echo "[Step 3] Intersecting with Genotypes (using IID)..."

python3 << EOF
import pandas as pd
import sys

triplets = pd.read_csv("$OUTDIR/tmp_triplets.csv")

# PLINK RAW: FID is 0, IID is the Subject ID (e.g. NEUCW292DYJ)
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])

# Expression: Normalize underscores to hyphens
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
expr.columns = [c.replace('_', '-') for c in expr.columns]

# PCA: Column 1 is IID
pca = pd.read_csv("$PCA", sep=r'\s+', header=None).iloc[:, [1,2,3,4,5,6]]
pca.columns = ['ID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']

# FINAL INTERSECTION
# Use 'subject' to match IID in genotypes and PCA
# Use 'hra' to match columns in Expression
overlap = triplets[
    triplets['subject'].isin(raw_df['IID']) & 
    triplets['hra'].isin(expr.columns) &
    triplets['subject'].isin(pca['ID'])
].copy()

print(f"  Final intersection count: {len(overlap)} samples.")
if overlap.empty:
    print("ERROR: Intersection is still empty. Check Subject ID strings.")
    sys.exit(1)

overlap.to_csv("$OUTDIR/tmp_final_overlap.csv", index=False)
EOF

##############################################
# STEP 4: Generate Expression Matrix
##############################################
echo "[Step 4] Finalizing expression file..."

python3 << EOF
import pandas as pd
overlap = pd.read_csv("$OUTDIR/tmp_final_overlap.csv")
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
expr.columns = [c.replace('_', '-') for c in expr.columns]

# Select samples and rename columns to DNA ID (hda) for Matrix eQTL
expr_out = expr[overlap['hra']]
expr_out.columns = overlap['hda']
expr_out.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')
EOF

##############################################
# STEP 5: Generate SNP Matrix
##############################################
echo "[Step 5] Finalizing SNP file..."

python3 << EOF
import pandas as pd
import numpy as np
overlap = pd.read_csv("$OUTDIR/tmp_final_overlap.csv")
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])

# Filter genotypes and set index to DNA ID (hda)
snp_out = raw_df[raw_df['IID'].isin(overlap['subject'])].copy()
# Map Subject back to HDA for the index
id_map = dict(zip(overlap['subject'], overlap['hda']))
snp_out['DNA_ID'] = snp_out['IID'].map(id_map)
snp_out = snp_out.set_index('DNA_ID').drop(columns=['FID', 'IID']).T

# MAF Filter
af = np.nanmean(snp_out.values, axis=1) / 2.0
maf = np.minimum(af, 1 - af)
keep_snps = np.where(np.nan_to_num(maf) >= 0.05)[0]
snp_out.iloc[keep_snps, :].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t')
EOF

##############################################
# STEP 6: Generate Covariate Matrix
##############################################
echo "[Step 6] Finalizing covariates..."

python3 << EOF
import pandas as pd
overlap = pd.read_csv("$OUTDIR/tmp_final_overlap.csv")
cov = pd.read_csv("$OUTDIR/tmp_clean_cov.tsv", sep='\t')
pca = pd.read_csv("$PCA", sep=r'\s+', header=None).iloc[:, [1,2,3,4,5,6]]
pca.columns = ['ID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']

cov_final = overlap.merge(cov, left_on='hra', right_on='externalsampleid')
cov_final = cov_final.merge(pca, left_on='subject', right_on='ID')

cov_final['sex_bin'] = cov_final['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})
cov_cols = ['sex_bin', 'age_at_death', 'rin', 'rna_skew', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
# Use hda as the column header for Matrix eQTL alignment
cov_final.set_index('hda')[cov_cols].T.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t')
EOF

##############################################
# STEP 7: Location Files
##############################################
echo "[Step 7] Generating location files..."
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt

echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3}' $GTF_BED6 >> $OUTDIR/gene_location.txt

##############################################
# STEP 8-10: Cleanup and Summary
##############################################
echo "[Step 8] Cleaning up..."
rm -f $OUTDIR/tmp_*

echo "============================================"
echo "  Prep Complete for $TISSUE"
echo "  $(date)"
echo "============================================"
