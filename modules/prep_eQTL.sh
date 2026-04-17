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
METADATA_CSV=/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv
BREAKDOWN_TSV=$QTL_DIR/patient_tissue_breakdown.tsv
COV=$QTL_DIR/covariates.tsv
RAW=$PLINK/joint_autosomes_matrixEQTL.raw
BIM=$PLINK/joint_autosomes_filtered_bed.bim
PCA=$PLINK/joint_pca.eigenvec
GTF_BED6=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6

echo "============================================"
echo "  eQTL Prep for $TISSUE"
echo "============================================"

python3 << EOF
import pandas as pd
import sys

# 1. LOAD DATA
breakdown = pd.read_csv("$BREAKDOWN_TSV", sep='\t')
meta = pd.read_csv("$METADATA_CSV", low_memory=False)

# CRITICAL FIX: Only load the necessary join key and unique columns from COV
# This prevents the name collision (externalsubjectid_x / _y)
cov_full = pd.read_csv("$COV", sep='\t')
cov_cols_to_keep = ['externalsampleid'] + [c for c in cov_full.columns if c not in meta.columns]
cov = cov_full[cov_cols_to_keep]

pca = pd.read_csv("$PCA", sep='\s+').rename(columns={'#IID': 'IID'})

# 2. FILTERING
patterns = {
    "Motor_Cortex": "Motor Cortex Lateral|Motor Cortex Medial|Lateral Motor Cortex|Medial Motor Cortex|Primary Motor Cortex L|Primary Motor Cortex M|Cortex_Motor_Unspecified|Cortex_Motor_BA4|BA4 Motor Cortex|Lateral_motor_cortex|Motor Cortex|BA4",
    "Cervical_Spinal_Cord": "Spinal_Cord_Cervical|Cervical Spinal Cord|Cervical_spinal_cord|Spinal_cord_Cervical|Cervical",
    "Lumbar_Spinal_Cord": "Lumbar Spinal Cord|Spinal_Cord_Lumbosacral|Lumbosacral_Spinal_Cord|Lumbar_spinal_cord|Lumbar|Lumbosacral",
    "Thoracic_Spinal_Cord": "Thoracic Spinal Cord|Thoracic",
    "Frontal_Cortex": "Frontal Cortex|Frontal",
    "Cerebellum": "Cerebellum"
}
pattern = patterns.get("$TISSUE", "$TISSUE")

approved_subjects = breakdown[breakdown['tissues'].str.contains("$TISSUE", na=False)]['subject_id'].unique()
meta_filt = meta[meta['externalsubjectid'].isin(approved_subjects)].copy()
mask = meta_filt.apply(lambda row: row.astype(str).str.contains(pattern, case=False).any(), axis=1)
meta_tissue = meta_filt[mask].copy()

# Select best replicate by RIN
def get_rin(row):
    try:
        val1 = float(row.iloc[16]) if pd.notnull(row.iloc[16]) else 0
        val2 = float(row.iloc[17]) if pd.notnull(row.iloc[17]) else 0
        return max(val1, val2)
    except: return 0

meta_tissue['RIN_score'] = meta_tissue.apply(get_rin, axis=1)
meta_unique = meta_tissue.sort_values('RIN_score', ascending=False).drop_duplicates(subset='externalsubjectid')

# 3. INTERSECTION
raw_df = pd.read_csv("$RAW", sep='\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
expr_headers = pd.read_csv("$EXPR", sep='\t', nrows=0).columns.tolist()

meta_unique['hra_clean'] = meta_unique['externalsampleid'].str.replace('-', '_')
common_subjects = sorted(list(set(meta_unique['externalsubjectid']) & set(raw_df['IID'])))

sub_to_hra = {}
final_aligned_subjects = []
for s in common_subjects:
    hra = meta_unique.loc[meta_unique['externalsubjectid'] == s, 'hra_clean'].values[0]
    matches = [c for c in expr_headers if c.startswith(hra)]
    if matches:
        final_aligned_subjects.append(s)
        sub_to_hra[s] = matches[0]

num_final = len(final_aligned_subjects)

# 4. ALIGN EXPRESSION & SNPS
expr = pd.read_csv("$EXPR", sep='\t', usecols=['gene_id'] + [sub_to_hra[s] for s in final_aligned_subjects], index_col=0)
expr_final = expr[[sub_to_hra[s] for s in final_aligned_subjects]]
expr_final.columns = final_aligned_subjects

snp_final = raw_df[raw_df['IID'].isin(final_aligned_subjects)].set_index('IID').drop(columns=['FID']).T
def clean_snp_id(full_id):
    parts = str(full_id).split(':')
    return f"chr{parts[0]}:{parts[1]}" if len(parts) >= 2 else full_id
snp_final.index = [clean_snp_id(idx) for idx in snp_final.index]

# 5. MERGE COVARIATES
meta_aligned = meta_unique[meta_unique['externalsubjectid'].isin(final_aligned_subjects)]

# First join: Metadata + Covariates (only externalsampleid is common)
cov_merged = pd.merge(meta_aligned, cov, on='externalsampleid')

# Second join: + PCA
cov_merged = pd.merge(cov_merged, pca, left_on='externalsubjectid', right_on='IID')

# Binary encoding
cov_merged['is_als'] = cov_merged['subject_group'].apply(lambda x: 1 if 'ALS' in str(x) else 0)
cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})

# Align order to expression matrix
cov_final = cov_merged.set_index('externalsubjectid').reindex(final_aligned_subjects)
cov_final = cov_final[['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']].T

# 6. SAVE
expr_final.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="geneid")
snp_final[final_aligned_subjects].to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="snpid")
cov_final.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', index=True, index_label="id", quoting=0)

print(f"SUCCESS: Processed {num_final} samples for $TISSUE")
EOF

### Location Files ###
echo "Generating location files..."
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt

echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3}' $GTF_BED6 >> $OUTDIR/gene_location.txt
