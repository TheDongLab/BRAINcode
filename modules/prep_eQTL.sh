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
echo "[Step 3] Intersecting with Genotypes and Expression..."

python3 << EOF
import pandas as pd
import sys

triplets = pd.read_csv("$OUTDIR/tmp_triplets.csv")

# 1. Load Genotypes
# We only care about the columns FID and IID here
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=[0, 1])

# 2. Load Expression & Normalize
# We will create a normalized list of headers for matching
expr_raw = pd.read_csv("$EXPR", sep='\t', index_col=0, nrows=1)
# Create a mapping of { 'CGND-HRA-00727' : 'CGND_HRA_00727' }
expr_mapping = {c.replace('_', '-'): c for c in expr_raw.columns}

# 3. INTERSECT
# Check if Subject is in PLINK IID AND normalized HRA is in Expression
overlap = triplets[
    triplets['subject'].isin(raw_df['IID']) & 
    triplets['hra'].isin(expr_mapping.keys())
].copy()

# 4. HANDLE PCA (The Broken File)
# Since your PCA file lacks IDs, we will assume it is ordered exactly 
# like your PLINK .raw file. If it has the same number of rows, we map them.
pca_raw = pd.read_csv("$PCA", sep=r'\s+', header=None)
if len(pca_raw) == len(raw_df):
    print("  PCA file matches PLINK row count. Mapping PCs by position.")
    # Map Subject IDs to PCs based on their position in the .raw file
    pca_cols = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']
    pca_data = pca_raw.iloc[:, 0:5] # Take first 5 columns
    pca_data.columns = pca_cols
    pca_data['ID'] = raw_df['IID'].values
    
    # Now check if our overlap samples have PCA data
    overlap = overlap[overlap['subject'].isin(pca_data['ID'])]
    pca_data.to_csv("$OUTDIR/tmp_fixed_pca.csv", index=False)
else:
    print(f"ERROR: PCA rows ({len(pca_raw)}) don't match PLINK rows ({len(raw_df)}).")
    sys.exit(1)

print(f"  Final intersection count: {len(overlap)} samples.")
if overlap.empty:
    print("ERROR: Intersection still 0. Check for hidden spaces in files.")
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
