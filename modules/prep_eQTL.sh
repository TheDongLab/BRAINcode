#!/bin/bash
#SBATCH --job-name=prep_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=649G

set -euo pipefail

module load BCFtools

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
COV_FILE=/home/zw529/donglab/data/target_ALS/QTL/covariates.tsv
SEX_MISMATCH_FILE=$OUTDIR/potential_sex_mismatches.txt

# --- GENOTYPE INPUTS ---
RAW=$PLINK/joint_all_chrs_matrixEQTL.raw
BIM=$PLINK/joint_all_chrs_filtered_bed.bim
PCA=$PLINK/joint_pca.eigenvec

# --- REFERENCE FILES ---
GTF_BED6=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6
MAP_FILE=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/ncbi_to_ucsc.txt
DBSNP_VCF=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/GCF_000001405.40.gz

echo "============================================"
echo "  eQTL Prep for $TISSUE"
echo "============================================"

# ─────────────────────────────────────────────────────────────────
# STEP A: HIGH-PERFORMANCE RSID EXTRACTION VIA BCFTOOLS (Path B)
# ─────────────────────────────────────────────────────────────────
echo "Streaming dbSNP variant IDs matching your dataset..."
# 1. Convert BIM coordinates to an internal 1-based region query file using the NC_ map
TMP_REGIONS=$OUTDIR/tmp_regions_ncbi.bed
awk -v map_file="$MAP_FILE" '
BEGIN {
    while ((getline < map_file) > 0) { map[$2] = $1 }
    close(map_file)
}
{
    ucsc = ($1 ~ /^chr/) ? $1 : "chr"$1;
    ncbi = (ucsc in map) ? map[ucsc] : ucsc;
    print ncbi"\t"$4"\t"$4
}' $BIM > $TMP_REGIONS

# 2. Extract matching variants from dbSNP using the Tabix index (instant, 0 RAM footprint)
TMP_RSID_MAP=$OUTDIR/tmp_coord_to_rsid.txt
bcftools query -R $TMP_REGIONS -f '%CHROM:%POS\t%ID\n' $DBSNP_VCF > $TMP_RSID_MAP

rm -f $TMP_REGIONS
echo "dbSNP query complete. Extracted subset mapping file."

# ─────────────────────────────────────────────────────────────────
# STEP B: PANDAS INTERSECTION & MATRIX GENERATION
# ─────────────────────────────────────────────────────────────────
python3 << EOF
import pandas as pd
import sys
import numpy as np

# 1. LOAD DATA & CHROMOSOME MAPS
breakdown = pd.read_csv("$BREAKDOWN_TSV", sep='\t')
meta = pd.read_csv("$METADATA_CSV", low_memory=False)
pca = pd.read_csv("$PCA", sep=r'\s+').rename(columns={'#IID': 'IID'})

# Load the reference dictionary to map UCSC 'chr' naming conventions to RefSeq 'NC_' tags
ucsc_to_ncbi = {}
with open("$MAP_FILE", 'r') as f:
    for line in f:
        if line.strip():
            ncbi, ucsc = line.strip().split('\t')
            ucsc_to_ncbi[ucsc] = ncbi

# Load the coordinate-to-rsID extracted mapping subset (Path B)
coord_to_rsid = {}
with open("$TMP_RSID_MAP", 'r') as f:
    for line in f:
        if line.strip():
            coord, rsid = line.strip().split('\t')
            # Only map if dbSNP actually returned a valid rsID string
            if rsid and rsid != '.':
                coord_to_rsid[coord] = rsid

# Load sex mismatch exclusions
try:
    with open("$SEX_MISMATCH_FILE") as f:
        sex_mismatch_samples = {
            line.strip() for line in f if line.strip()
        }
except FileNotFoundError:
    sex_mismatch_samples = set()

print(f"DEBUG: Loaded {len(sex_mismatch_samples)} sex-mismatch samples")

# Remove from metadata immediately
meta = meta[
    ~meta['externalsampleid'].isin(sex_mismatch_samples)
].copy()

# 2. DEFINE TISSUE PATTERN
patterns = {
    "Motor_Cortex": "Motor Cortex Lateral|Motor Cortex Medial|Lateral Motor Cortex|Medial Motor Cortex|Primary Motor Cortex L|Primary Motor Cortex M|Cortex_Motor_Unspecified|Cortex_Motor_BA4|BA4 Motor Cortex|Lateral_motor_cortex|Motor Cortex|BA4",
    "Cervical_Spinal_Cord": "Spinal_Cord_Cervical|Cervical Spinal Cord|Cervical_spinal_cord|Spinal_cord_Cervical|Cervical",
    "Lumbar_Spinal_Cord": "Lumbar Spinal Cord|Spinal_Cord_Lumbosacral|Lumbosacral_Spinal_Cord|Lumbar_spinal_cord|Lumbar|Lumbosacral",
    "Thoracic_Spinal_Cord": "Thoracic Spinal Cord|Thoracic",
    "Frontal_Cortex": "Frontal Cortex|Frontal",
    "Cerebellum": "Cerebellum"
}
pattern = patterns.get("$TISSUE", "$TISSUE")

# 3. FILTERING & BEST REPLICATE SELECTION
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

# FILTER: PMI <= 40 hours & RIN >= 3
print(f"DEBUG: Before filtering: {len(meta_tissue)} samples")
meta_tissue = meta_tissue[(meta_tissue['post_mortem_interval_in_hours'] <= 40) & (meta_tissue['RIN_score'] >= 3)]
print(f"DEBUG: After PMI/RIN filter: {len(meta_tissue)} samples")

meta_unique = meta_tissue.sort_values('RIN_score', ascending=False).drop_duplicates(subset='externalsubjectid')

# 4. TRIPLE INTERSECTION
raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
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

# 5. ALIGN EXPRESSION & SNPS
expr = pd.read_csv("$EXPR", sep='\t', usecols=['gene_id'] + [sub_to_hra[s] for s in final_aligned_subjects], index_col=0)
expr_final = expr[[sub_to_hra[s] for s in final_aligned_subjects]]
expr_final.columns = final_aligned_subjects

# Correct handling of SNP indices to preserve additive encoding (0,1,2)
snp_matrix = raw_df[raw_df['IID'].isin(final_aligned_subjects)].set_index('IID').drop(columns=['FID'])
snp_matrix = snp_matrix.loc[final_aligned_subjects]
snp_final = snp_matrix.T

def convert_to_rsid(full_id):
    # Strip recode suffix (_A, _G)
    clean_id = full_id.rsplit('_', 1)[0] if '_' in full_id else full_id
    parts = str(clean_id).split(':')
    if len(parts) >= 2:
        chrom, pos = parts[0], parts[1]
        ucsc_chrom = chrom if chrom.startswith('chr') else f"chr{chrom}"
        ncbi_chrom = ucsc_to_ncbi.get(ucsc_chrom, ucsc_chrom)
        
        # Look up coordinate string in our dynamic dbSNP dictionary
        coord_key = f"{ncbi_chrom}:{pos}"
        return coord_to_rsid.get(coord_key, coord_key) # Fallback to coordinate if variant isn't in dbSNP
    return clean_id

snp_final.index = [convert_to_rsid(idx) for idx in snp_final.index]

# CRITICAL FIX: Remove duplicate SNP IDs (keep first occurrence only)
snp_final = snp_final[~snp_final.index.duplicated(keep='first')]

print(f"DEBUG: SNP matrix shape after rsID conversion: {snp_final.shape}")

# 6. COVARIATE MERGE
cov_full = pd.read_csv("$COV_FILE", sep='\t')
cov_cols_to_keep = ['externalsampleid'] + [c for c in cov_full.columns if c not in meta.columns]
cov = cov_full[cov_cols_to_keep].drop_duplicates(subset='externalsampleid')

meta_aligned = meta_unique[meta_unique['externalsubjectid'].isin(final_aligned_subjects)].copy()
cov_merged = pd.merge(meta_aligned, cov, on='externalsampleid', how='inner')
cov_merged = pd.merge(cov_merged, pca, left_on='externalsubjectid', right_on='IID', how='inner')
cov_merged = cov_merged.sort_values('RIN_score', ascending=False).drop_duplicates(subset='externalsubjectid')

cov_merged['is_als'] = cov_merged['subject_group'].apply(lambda x: 1 if 'ALS' in str(x) else 0)
cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})

cov_final = cov_merged.set_index('externalsubjectid').reindex(final_aligned_subjects)
cov_final = cov_final[['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']].T

# 7. SAVE
expr_final.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="geneid")
snp_final.to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="snpid")
cov_final.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', index=True, index_label="id", quoting=0)

print(f"SUCCESS: Processed {num_final} unique samples for $TISSUE")
EOF

# ─────────────────────────────────────────────────────────────────
# STEP C: GENERATE ALIGNED LOCATION FILES
# ─────────────────────────────────────────────────────────────────
echo "Generating deduplicated location files for $TISSUE..."

# 1. SNP Location: Converts chromosome format AND uses rsIDs matching snp_${TISSUE_DIR}.txt
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk -v ncbi_map="$MAP_FILE" -v rsid_map="$TMP_RSID_MAP" '
BEGIN {
    OFS="\t"
    # Load UCSC to NCBI map
    while ((getline < ncbi_map) > 0) { n_map[$2] = $1 }
    close(ncbi_map)
    # Load NCBI:POS to rsID map
    while ((getline < rsid_map) > 0) { r_map[$1] = $2 }
    close(rsid_map)
}
{
    ucsc = ($1 ~ /^chr/) ? $1 : "chr"$1;
    ncbi = (ucsc in n_map) ? n_map[ucsc] : ucsc;
    
    coord_key = ncbi":"$4;
    final_id = (coord_key in r_map) ? r_map[coord_key] : coord_key;
    
    print final_id, ncbi, $4;
}' $BIM | sort -u -k1,1 >> $OUTDIR/snp_location.txt

# 2. Gene Location: Extract from GTF, then force unique IDs based on column 1
echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3}' $GTF_BED6 | sort -u -k1,1 >> $OUTDIR/gene_location.txt

# Post-execution cleanup of local session mapping file
rm -f $TMP_RSID_MAP
echo "Location files complete for $TISSUE."
