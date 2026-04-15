#!/bin/bash
#SBATCH --job-name=prep_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=399G

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
# STEP 1: Pre-process Covariates
##############################################
echo "[Step 1] Cleaning covariates and encoding phenotype..."

python3 << EOF
import pandas as pd
df = pd.read_csv("$COV", sep='\t')
# Filter by quality metrics
bad_samples = df[df['rin'].isna() & (df['rna_skew'] > 1.0)]['externalsampleid'].tolist()
df = df[~df['externalsampleid'].isin(bad_samples)]

# Binary Encoding for Case/Control
def collapse(g):
    if 'ALS' in str(g): return 1
    if 'Non-Neurological' in str(g): return 0
    return None

df['is_als'] = df['subject_group'].apply(collapse)
df = df[df['is_als'].notna()]
df.to_csv("$OUTDIR/tmp_clean_cov.tsv", sep='\t', index=False)
EOF

##############################################
# STEP 2: ID Alignment (The "Bridge" Step)
##############################################
echo "[Step 2] Bridging DNA and RNA IDs in Python..."

python3 << EOF
import pandas as pd
import numpy as np
import sys

# Load Sample Map (Tabs verified: HDA_ID, Subject_ID)
map_df = pd.read_csv("$SAMPLE_MAP", sep='\t', header=None, names=['hda', 'subject'])
map_df['hda'] = map_df['hda'].astype(str).str.replace(r'-b\d+$', '', regex=True).str.strip()
map_df['subject'] = map_df['subject'].astype(str).str.strip()

# Load Metadata (Commas verified: RNA_ID, Subject_ID, Tissue)
# We read it without a header first to handle the specific structure you provided
meta_df = pd.read_csv("$META", header=None, names=['hra', 'subject', 'tissue'])
meta_filt = meta_df[meta_df['tissue'] == "$TISSUE"].copy()
meta_filt['hra'] = meta_filt['hra'].astype(str).str.strip()
meta_filt['subject'] = meta_filt['subject'].astype(str).str.strip()

# Create the Triplets Bridge
triplets = pd.merge(meta_filt, map_df, on='subject')
triplets.to_csv("$OUTDIR/tmp_triplets.csv", index=False)

print(f"  Mapped {len(triplets)} candidates for $TISSUE.")
if triplets.empty:
    print("ERROR: Step 2 failed. No overlap between metadata subjects and sample map.")
    sys.exit(1)
EOF

##############################################
# STEP 3: Load and Intersect High-D Data
##############################################
echo "[Step 3] Intersecting candidates with RAW, EXPR, and PCA..."

python3 << EOF
import pandas as pd
import numpy as np
import sys

triplets = pd.read_csv("$OUTDIR/tmp_triplets.csv")

# Load RAW genotypes (Use raw string for sep)
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
raw_df['FID'] = raw_df['FID'].astype(str).str.replace(r'-b\d+$', '', regex=True)

# Load Expression and normalize headers (Underscore -> Hyphen)
expr = pd.read_csv("$EXPR", sep='\t', index_col=0)
expr.columns = [c.replace('_', '-') for c in expr.columns]

# Load PCA
pca = pd.read_csv("$PCA", sep=r'\s+', header=None).iloc[:, [1,2,3,4,5,6]]
pca.columns = ['ID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
pca['ID'] = pca['ID'].astype(str).str.replace(r'-b\d+$', '', regex=True)

# Final Intersection
overlap = triplets[
    triplets['hda'].isin(raw_df['FID']) & 
    triplets['hra'].isin(expr.columns) &
    triplets['hda'].isin(pca['ID'])
].copy()

print(f"  Final intersection count: {len(overlap)}")
if overlap.empty:
    print("ERROR: Intersection is empty. DNA IDs likely mismatch between SampleMap and PLINK.")
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

expr_out = expr[overlap['hra']]
expr_out.columns = overlap['hda'] # Column names MUST be DNA IDs for Matrix eQTL
expr_out.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')
EOF

##############################################
# STEP 5: Generate SNP Matrix
##############################################
echo "[Step 5] Finalizing SNP file (MAF > 0.05)..."

python3 << EOF
import pandas as pd
import numpy as np
overlap = pd.read_csv("$OUTDIR/tmp_final_overlap.csv")
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
raw_df['FID'] = raw_df['FID'].astype(str).str.replace(r'-b\d+$', '', regex=True)

snp_out = raw_df[raw_df['FID'].isin(overlap['hda'])].set_index('FID').drop(columns=['IID']).T
# Filter for Minor Allele Frequency
af = np.nanmean(snp_out.values, axis=1) / 2.0
maf = np.minimum(af, 1 - af)
keep_snps = np.where(np.nan_to_num(maf) >= 0.05)[0]
snp_out.iloc[keep_snps, :].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t')
EOF

##############################################
# STEP 6: Generate Covariate Matrix
##############################################
echo "[Step 6] Finalizing covariate and PCA file..."

python3 << EOF
import pandas as pd
overlap = pd.read_csv("$OUTDIR/tmp_final_overlap.csv")
cov = pd.read_csv("$OUTDIR/tmp_clean_cov.tsv", sep='\t')
pca = pd.read_csv("$PCA", sep=r'\s+', header=None).iloc[:, [1,2,3,4,5,6]]
pca.columns = ['ID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
pca['ID'] = pca['ID'].astype(str).str.replace(r'-b\d+$', '', regex=True)

# Merge metadata with original covariates and PCA
cov_final = overlap.merge(cov, left_on='hra', right_on='externalsampleid')
cov_final = cov_final.merge(pca, left_on='hda', right_on='ID')

cov_final['sex_bin'] = cov_final['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})
# Transpose for Matrix eQTL format
cov_cols = ['sex_bin', 'age_at_death', 'rin', 'rna_skew', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
cov_final.set_index('hda')[cov_cols].T.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t')
EOF

##############################################
# STEP 7: Generate Location Files
##############################################
echo "[Step 7] Generating SNP and Gene location files..."

# SNP Loc
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt

# Gene Loc
echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {
    gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3
}' $GTF_BED6 >> $OUTDIR/gene_location.txt

##############################################
# STEP 8: Cleanup and Validation
##############################################
echo "[Step 8] Cleaning up and validating output..."

rm -f $OUTDIR/tmp_*

python3 << EOF
import pandas as pd
# Check that all 3 main files have the same number of columns (samples)
expr = pd.read_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t', nrows=1)
snp = pd.read_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t', nrows=1)
cov = pd.read_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', nrows=1)

print(f"  Samples in Expression: {len(expr.columns)-1}")
print(f"  Samples in SNPs      : {len(snp.columns)-1}")
print(f"  Samples in Covariates: {len(cov.columns)-1}")

if not (len(expr.columns) == len(snp.columns) == len(cov.columns)):
    print("WARNING: Sample counts do not match! Check column alignment.")
EOF

echo "============================================"
echo "  Success: $TISSUE Prep Complete"
echo "  Location: $OUTDIR"
echo "  $(date)"
echo "============================================"
