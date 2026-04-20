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
PLINK_PREFIX=$QTL_DIR/plink/joint_autosomes_filtered_bed
EXPR=$QTL_DIR/expression_matrix.txt
METADATA_CSV=/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv
BREAKDOWN_TSV=$QTL_DIR/patient_tissue_breakdown.tsv
COV_FILE=/home/zw529/donglab/data/target_ALS/QTL/covariates.tsv
BIM=${PLINK_PREFIX}.bim
PCA=$QTL_DIR/plink/joint_pca.eigenvec
GTF_BED6=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6

# Define intermediate raw file for this run
RAW_GENOTYPES=$OUTDIR/genotypes_additive.raw

echo "============================================"
echo "  eQTL Prep for $TISSUE"
echo "============================================"

# STEP 0: Generate fresh Additive Genotypes to ensure 0, 1, 2 encoding
echo "Exporting additive genotypes from PLINK..."
plink2 --bfile $PLINK_PREFIX \
       --recode A \
       --out ${RAW_GENOTYPES%.raw}

python3 << EOF
import pandas as pd
import numpy as np

# 1. LOAD DATA
meta = pd.read_csv("$METADATA_CSV", low_memory=False)
breakdown = pd.read_csv("$BREAKDOWN_TSV", sep='\t')
pca = pd.read_csv("$PCA", sep='\s+').rename(columns={'#IID': 'IID'})

# 2. DEFINE TISSUE PATTERN
patterns = {
    "Motor_Cortex": "Motor Cortex Lateral|Motor Cortex Medial|Lateral Motor Cortex|Medial Motor Cortex|Primary Motor Cortex L|Primary Motor Cortex M|Cortex_Motor_Unspecified|Cortex_Motor_BA4|BA4 Motor Cortex|Lateral_motor_cortex|Motor Cortex|BA4",
    "Cerebellum": "Cerebellum"
}
pattern = patterns.get("$TISSUE", "$TISSUE")

# 3. FILTERING & SUBJECT ALIGNMENT
approved_subjects = breakdown[breakdown['tissues'].str.contains("$TISSUE", na=False)]['subject_id'].unique()
meta_filt = meta[meta['externalsubjectid'].isin(approved_subjects)].copy()
mask = meta_filt.apply(lambda row: row.astype(str).str.contains(pattern, case=False).any(), axis=1)
meta_tissue = meta_filt[mask].copy()

def get_rin(row):
    try:
        val1 = float(row.iloc[16]) if pd.notnull(row.iloc[16]) else 0
        val2 = float(row.iloc[17]) if pd.notnull(row.iloc[17]) else 0
        return max(val1, val2)
    except: return 0

meta_tissue['RIN_score'] = meta_tissue.apply(get_rin, axis=1)
meta_unique = meta_tissue.sort_values('RIN_score', ascending=False).drop_duplicates(subset='externalsubjectid')

# 4. LOAD GENOTYPES & INTERSECT
raw_df = pd.read_csv("$RAW_GENOTYPES", sep='\s+')
common_subjects = sorted(list(set(meta_unique['externalsubjectid']) & set(raw_df['IID'])))

# ALIGN SAMPLES
expr_headers = pd.read_csv("$EXPR", sep='\t', nrows=0).columns.tolist()
meta_unique['hra_clean'] = meta_unique['externalsampleid'].str.replace('-', '_')

final_subjects = []
sub_to_hra = {}
for s in common_subjects:
    hra = meta_unique.loc[meta_unique['externalsubjectid'] == s, 'hra_clean'].values[0]
    matches = [c for c in expr_headers if c.startswith(hra)]
    if matches:
        final_subjects.append(s)
        sub_to_hra[s] = matches[0]

# 5. PROCESS GENOTYPES (The "2" preserver)
# Subset to aligned subjects and drop PLINK metadata
snp_matrix = raw_df[raw_df['IID'].isin(final_subjects)].set_index('IID').drop(columns=['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'])

# Transpose
snp_final = snp_matrix.T

# Clean SNP IDs (Remove the _A suffix PLINK adds)
# e.g., "chr1:12345_A" -> "chr1:12345"
snp_final.index = [i.rsplit('_', 1)[0] if '_' in str(i) else i for i in snp_final.index]

# DEBUG: Verify we have 2s
val_counts = np.unique(snp_final.values[~np.isnan(snp_final.values)])
print(f"DEBUG: Found genotype values: {val_counts}")

# 6. PROCESS EXPRESSION
expr = pd.read_csv("$EXPR", sep='\t', usecols=['gene_id'] + [sub_to_hra[s] for s in final_subjects], index_col=0)
expr_final = expr[[sub_to_hra[s] for s in final_subjects]]
expr_final.columns = final_subjects

# 7. PROCESS COVARIATES
cov_full = pd.read_csv("$COV_FILE", sep='\t')
cov_cols = ['externalsampleid'] + [c for c in cov_full.columns if c not in meta.columns]
cov = cov_full[cov_cols].drop_duplicates(subset='externalsampleid')

meta_aligned = meta_unique[meta_unique['externalsubjectid'].isin(final_subjects)]
cov_merged = pd.merge(meta_aligned, cov, on='externalsampleid', how='inner')
cov_merged = pd.merge(cov_merged, pca, left_on='externalsubjectid', right_on='IID', how='inner')

cov_merged['is_als'] = cov_merged['subject_group'].apply(lambda x: 1 if 'ALS' in str(x) else 0)
cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})

cov_final = cov_merged.set_index('externalsubjectid').reindex(final_subjects)
cov_final = cov_final[['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']].T

# 8. SAVE (Force tab-separation and no quoting to keep files clean)
expr_final.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')
snp_final[final_subjects].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t')
cov_final.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t')

print(f"SUCCESS: {len(final_subjects)} samples processed.")
EOF

### Location Files ###
echo "Generating location files..."
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt
echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3}' $GTF_BED6 >> $OUTDIR/gene_location.txt
