#!/bin/bash
#SBATCH --job-name=prep_sQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.err
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=699G

set -euo pipefail

TISSUE="$1" 
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')
OUTDIR=/home/zw529/donglab/data/target_ALS/$TISSUE_DIR/sQTL
mkdir -p $OUTDIR

### Paths ###
QTL_DIR=/home/zw529/donglab/data/target_ALS/QTL
PLINK=$QTL_DIR/plink
SPLICING_RAW=$QTL_DIR/splicing_matrix.txt 
SPLICING_LOC_SRC=$OUTDIR/tmp_splicing_events_hg38.bed

METADATA_CSV=/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv
BREAKDOWN_TSV=$QTL_DIR/patient_tissue_breakdown.tsv
COV_FILE=/home/zw529/donglab/data/target_ALS/QTL/covariates.tsv
SEX_MISMATCH_FILE=/home/zw529/donglab/data/target_ALS/$TISSUE_DIR/eQTL/potential_sex_mismatches.txt

### ACTIVE ALL-CHROMOSOME IMPUTED INPUTS ###
RAW=$PLINK/joint_all_chrs_matrixEQTL.raw
BIM=$PLINK/joint_all_chrs_filtered_bed.bim
PCA=$PLINK/joint_pca.eigenvec

echo "============================================"
echo "  sQTL Prep for $TISSUE"
echo "============================================"

##############################################
# STEP 0: Generate BED from Matrix IDs
##############################################
echo "[0] Generating BED file from matrix IDs..."
python3 << EOF
import sys

try:
    encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    opened = False
    
    for enc in encodings:
        try:
            with open("$SPLICING_RAW", 'r', encoding=enc, errors='replace') as f:
                header = f.readline()  # Skip header
                with open("$SPLICING_LOC_SRC", 'w') as out:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        parts = line.split('\t')
                        if len(parts) < 1:
                            continue
                        event_id = parts[0]
                        event_id = ''.join(c for c in event_id if ord(c) < 128)
                        
                        try:
                            if ':' in event_id and '-' in event_id:
                                chr_part, pos_part = event_id.rsplit(':', 1)
                                start, end = pos_part.split('-')
                                out.write(f"{chr_part}\t{start}\t{end}\t{event_id}\n")
                        except (ValueError, IndexError):
                            pass
            opened = True
            print(f"Successfully read splicing matrix with encoding: {enc}", file=sys.stderr)
            break
        except UnicodeDecodeError:
            continue
    
    if not opened:
        raise ValueError("Could not read splicing matrix with any encoding")
except Exception as e:
    print(f"ERROR in STEP 0: {e}", file=sys.stderr)
    sys.exit(1)

print("STEP 0: BED file generated successfully", file=sys.stderr)
EOF

##############################################
# STEP 1-6: Load, Filter, and Align Matrix Data
##############################################
echo "[1] Processing cohort alignments..."
python3 << EOF
import pandas as pd
import sys
import numpy as np

try:
    # 1. LOAD DATA
    breakdown = pd.read_csv("$BREAKDOWN_TSV", sep='\t')
    meta = pd.read_csv("$METADATA_CSV", low_memory=False)
    pca = pd.read_csv("$PCA", sep=r'\s+').rename(columns={'#IID': 'IID'})

    # Load sex mismatch exclusions
    try:
        with open("$SEX_MISMATCH_FILE") as f:
            sex_mismatch_samples = {line.strip() for line in f if line.strip()}
    except FileNotFoundError:
        sex_mismatch_samples = set()

    print(f"DEBUG: Loaded {len(sex_mismatch_samples)} sex-mismatch samples")

    # Remove from metadata immediately
    meta = meta[~meta['externalsampleid'].isin(sex_mismatch_samples)].copy()

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
    
    encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    splicing_headers = []
    for enc in encodings:
        try:
            splicing_headers = pd.read_csv("$SPLICING_RAW", sep='\t', nrows=0, encoding=enc).columns.tolist()
            break
        except UnicodeDecodeError:
            continue

    meta_unique['hra_clean'] = meta_unique['externalsampleid'].str.replace('-', '_')
    common_subjects = sorted(list(set(meta_unique['externalsubjectid']) & set(raw_df['IID'])))

    sub_to_hra = {}
    final_aligned_subjects = []
    for s in common_subjects:
        hra = meta_unique.loc[meta_unique['externalsubjectid'] == s, 'hra_clean'].values[0]
        matches = [c for c in splicing_headers if c.startswith(hra)]
        if matches:
            final_aligned_subjects.append(s)
            sub_to_hra[s] = matches[0]

    num_final = len(final_aligned_subjects)

    # 5. ALIGN SPLICING & SNPS
    splicing_df = None
    for enc in encodings:
        try:
            splicing_df = pd.read_csv("$SPLICING_RAW", sep='\t', usecols=['chrom'] + [sub_to_hra[s] for s in final_aligned_subjects], index_col=0, encoding=enc)
            break
        except UnicodeDecodeError:
            continue

    splicing_final = splicing_df[[sub_to_hra[s] for s in final_aligned_subjects]]
    splicing_final.columns = final_aligned_subjects

    # Non-Zero Variance Row Filter
    row_vars = splicing_final.var(axis=1, skipna=True)
    splicing_final = splicing_final[row_vars > 1e-10]

    # Correct handling of SNP indices to preserve additive encoding (0,1,2)
    snp_matrix = raw_df[raw_df['IID'].isin(final_aligned_subjects)].set_index('IID').drop(columns=['FID'])
    snp_matrix = snp_matrix.loc[final_aligned_subjects]
    snp_final = snp_matrix.T

    def clean_snp_id(full_id):
        clean_id = full_id.rsplit('_', 1)[0] if '_' in full_id else full_id
        parts = str(clean_id).split(':')
        if len(parts) >= 2:
            return clean_id if str(parts[0]).startswith('chr') else f"chr{parts[0]}:{parts[1]}"
        return clean_id

    snp_final.index = [clean_snp_id(idx) for idx in snp_final.index]
    snp_final = snp_final[~snp_final.index.duplicated(keep='first')]

    print(f"DEBUG: SNP matrix shape after dedup: {snp_final.shape}")
    print(f"DEBUG: Sample order check - first 3 samples: {list(snp_final.columns[:3])}")

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
    splicing_final.to_csv("$OUTDIR/splicing_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="geneid")
    snp_final.to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="snpid")
    cov_final.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', index=True, index_label="id", quoting=0)

    print(f"SUCCESS: Processed {num_final} unique samples for $TISSUE")

except Exception as e:
    print(f"ERROR in processing steps: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF

##############################################
# STEP 7: Coordinate Mapping Location Files
##############################################
echo "Generating deduplicated location files for $TISSUE..."

# 1. SNP Location file generation
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {
    chrom = ($1 ~ /^chr/) ? $1 : "chr"$1;
    print "chr"$1":"$4, chrom, $4;
}' $BIM | sed 's/chrchr/chr/g' | sort -u -k1,1 >> $OUTDIR/snp_location.txt

# 2. Splicing target coordinate maps
echo -e "geneid\tchr\tleft\tright" > $OUTDIR/splicing_location.txt
awk 'BEGIN{OFS="\t"} { print $4, $1, $2, $3 }' "$SPLICING_LOC_SRC" | sort -u -k1,1 >> $OUTDIR/splicing_location.txt

# Clean up temporary structural coordinates
rm -f "$SPLICING_LOC_SRC"

echo "Location files complete for $TISSUE."

##############################################
# STEP 8: Final Summary Output
##############################################
echo "============================================"
echo "  sQTL Prep complete for $TISSUE"
echo "  Output directory: $OUTDIR"
echo "  $(date)"
echo "============================================"

# List output files to mirror initial script diagnostics
echo "Output files:"
ls -lh "$OUTDIR"/*.txt 2>/dev/null || echo "No output files found"
