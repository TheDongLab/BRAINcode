#!/bin/bash
#SBATCH --job-name=prep_sQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.err
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=649G

set -euo pipefail

module load BCFtools

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

# --- GENOTYPE INPUTS ---
RAW=$PLINK/joint_all_chrs_matrixEQTL.raw
BIM=$PLINK/joint_all_chrs_filtered_bed.bim
PCA=$PLINK/joint_pca.eigenvec

# --- REFERENCE FILES ---
MAP_FILE=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/ncbi_to_ucsc.txt
DBSNP_VCF=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/GCF_000001405.40.gz

echo "============================================"
echo "  sQTL Prep for $TISSUE"
echo "============================================"

# ─────────────────────────────────────────────────────────────────
# STEP A: HIGH-PERFORMANCE RSID EXTRACTION VIA BCFTOOLS -T (Path B)
# ─────────────────────────────────────────────────────────────────
echo "Generating target position filter list..."
TMP_TARGETS=$OUTDIR/tmp_targets_ncbi.txt

# 1. Convert BIM positions into a clean, 3-column tab-delimited filter target file
awk -v map_file="$MAP_FILE" '
BEGIN {
    OFS="\t"
    while ((getline < map_file) > 0) { map[$2] = $1 }
    close(map_file)
}
{
    ucsc = ($1 ~ /^chr/) ? $1 : "chr"$1;
    ncbi = (ucsc in map) ? map[ucsc] : ucsc;
    print ncbi, $4, $4
}' $BIM > $TMP_TARGETS

echo "Streaming dbSNP variant IDs matching your dataset via row filtering..."
TMP_RSID_MAP=$OUTDIR/tmp_coord_to_rsid.txt

# 2. Extract matching variants from dbSNP using the non-hanging sequential text filter (-T)
bcftools query -T $TMP_TARGETS -f '%CHROM:%POS\t%ID\n' $DBSNP_VCF > $TMP_RSID_MAP

# Clean up the intermediate targets file
rm -f $TMP_TARGETS
echo "dbSNP query complete. Extracted subset mapping file successfully."

# ─────────────────────────────────────────────────────────────────
# STEP 0: Generate BED from Matrix IDs
# ─────────────────────────────────────────────────────────────────
echo "[0] Generating BED file from matrix IDs..."
python3 << EOF
import sys
try:
    for enc in ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']:
        try:
            with open("$SPLICING_RAW", 'r', encoding=enc, errors='replace') as f, open("$SPLICING_LOC_SRC", 'w') as out:
                f.readline()
                for line in f:
                    if not line.strip(): continue
                    event_id = ''.join(c for c in line.split('\t')[0] if ord(c) < 128)
                    if ':' in event_id and '-' in event_id:
                        chr_part, pos_part = event_id.rsplit(':', 1)
                        start, end = pos_part.split('-')
                        out.write(f"{chr_part.split(':')[0]}\t{start}\t{end}\t{event_id}\n")
            print(f"STEP 0: Complete using encoding: {enc}", file=sys.stderr)
            break
        except UnicodeDecodeError: continue
except Exception as e:
    print(f"ERROR in STEP 0: {e}", file=sys.stderr); sys.exit(1)
EOF

# ─────────────────────────────────────────────────────────────────
# STEP 1-6: Load, Filter, and Align Matrix Data
# ─────────────────────────────────────────────────────────────────
echo "[1] Processing cohort alignments with on-the-fly read filtering..."
python3 << EOF
import pandas as pd
import sys, os
import numpy as np
from collections import Counter

try:
    # 1. LOAD DATA & REMOVE SEX MISMATCHES
    breakdown = pd.read_csv("$BREAKDOWN_TSV", sep='\t')
    meta = pd.read_csv("$METADATA_CSV", low_memory=False)
    pca = pd.read_csv("$PCA", sep=r'\s+').rename(columns={'#IID': 'IID'})

    ucsc_to_ncbi = {}
    with open("$MAP_FILE", 'r') as f:
        for line in f:
            if line.strip():
                ncbi, ucsc = line.strip().split('\t')
                ucsc_to_ncbi[ucsc] = ncbi

    coord_to_rsid = {}
    with open("$TMP_RSID_MAP", 'r') as f:
        for line in f:
            if line.strip():
                coord, rsid = line.strip().split('\t')
                if rsid and rsid != '.':
                    coord_to_rsid[coord] = rsid

    try:
        with open("$SEX_MISMATCH_FILE") as f: sex_mismatch = {l.strip() for l in f if l.strip()}
    except FileNotFoundError: sex_mismatch = set()
    meta = meta[~meta['externalsampleid'].isin(sex_mismatch)].copy()

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

    # 3. FILTERING & BEST REPLICATE SELECTION (PMI <= 40 & RIN >= 3)
    approved_subjects = breakdown[breakdown['tissues'].str.contains("$TISSUE", na=False)]['subject_id'].unique()
    meta_filt = meta[meta['externalsubjectid'].isin(approved_subjects)].copy()
    meta_tissue = meta_filt[meta_filt.apply(lambda r: r.astype(str).str.contains(pattern, case=False).any(), axis=1)].copy()

    def get_rin(r):
        try: return max(float(r.iloc[16]) if pd.notnull(r.iloc[16]) else 0, float(r.iloc[17]) if pd.notnull(r.iloc[17]) else 0)
        except: return 0
    meta_tissue['RIN_score'] = meta_tissue.apply(get_rin, axis=1)
    meta_tissue = meta_tissue[(meta_tissue['post_mortem_interval_in_hours'] <= 40) & (meta_tissue['RIN_score'] >= 3)]
    meta_unique = meta_tissue.sort_values('RIN_score', ascending=False).drop_duplicates(subset='externalsubjectid').copy()

    # 4. TRIPLE INTERSECTION & INDIVIDUAL FILE READ-COUNT SCAN
    raw_df = pd.read_csv("$RAW", sep=r'\s+', usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE'])
    encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    splicing_headers = []
    for enc in encodings:
        try: splicing_headers = pd.read_csv("$SPLICING_RAW", sep='\t', nrows=0, encoding=enc).columns.tolist(); break
        except UnicodeDecodeError: continue

    meta_unique['hra_clean'] = meta_unique['externalsampleid'].str.replace('-', '_')
    common_subjects = sorted(list(set(meta_unique['externalsubjectid']) & set(raw_df['IID'])))

    sub_to_hra, final_aligned_subjects = {}, []
    junc_counts = Counter()
    rep_dir_base = "/home/zw529/donglab/data/target_ALS/" + "$TISSUE_DIR" + "/RNAseq/Processed"

    for s in common_subjects:
        row = meta_unique[meta_unique['externalsubjectid'] == s]
        hra = row['hra_clean'].values[0]
        raw_sample = row['externalsampleid'].values[0] 
        raw_sample_clean = raw_sample.replace('-', '_')
        
        matches = [c for c in splicing_headers if c.startswith(hra)]
        if matches:
            tsv_path = os.path.join(rep_dir_base, raw_sample_clean, "leafcutter", "psi", f"{raw_sample_clean}.leafcutter.PSI.tsv")
            if os.path.exists(tsv_path):
                final_aligned_subjects.append(s)
                sub_to_hra[s] = matches[0]
                try:
                    df_rep = pd.read_csv(tsv_path, sep=r'\s+', usecols=['chrom', 'strand', 'start', 'end', 'reads'])
                    df_rep['j_id'] = df_rep['chrom'].astype(str) + ':' + df_rep['strand'].astype(str) + ':' + df_rep['start'].astype(str) + '-' + df_rep['end'].astype(str)
                    junc_counts.update(df_rep.loc[df_rep['reads'] >= 5, 'j_id'].unique())
                except: pass

    if len(final_aligned_subjects) == 0:
        raise ValueError(f"CRITICAL ERROR: Zero sample files matched! Checked base: {rep_dir_base}")

    min_pass = max(1, int(len(final_aligned_subjects) * 0.10))
    valid_junctions = {j for j, c in junc_counts.items() if c >= min_pass}
    print(f"DEBUG: Mapped {len(final_aligned_subjects)} samples. Retained {len(valid_junctions)} junctions matching expression filters.")

    # 5. ALIGN SPLICING & SNPS
    splicing_df = None
    for enc in encodings:
        try:
            target_cols = ['junction_id'] + [sub_to_hra[s] for s in final_aligned_subjects]
            splicing_df = pd.read_csv("$SPLICING_RAW", sep='\t', usecols=target_cols, index_col=0, encoding=enc)
            break
        except UnicodeDecodeError: continue

    splicing_final = splicing_df[[sub_to_hra[s] for s in final_aligned_subjects]]
    splicing_final.columns = final_aligned_subjects
    splicing_final = splicing_final.loc[splicing_final.index.isin(valid_junctions)]
    splicing_final = splicing_final[splicing_final.var(axis=1, skipna=True) > 1e-10]

    snp_matrix = raw_df[raw_df['IID'].isin(final_aligned_subjects)].set_index('IID').drop(columns=['FID']).loc[final_aligned_subjects]
    snp_final = snp_matrix.T

    # Map SNP coordinate IDs dynamically to rsIDs
    def convert_to_rsid(full_id):
        clean_id = full_id.rsplit('_', 1)[0] if '_' in full_id else full_id
        parts = str(clean_id).split(':')
        if len(parts) >= 2:
            chrom, pos = parts[0], parts[1]
            ucsc_chrom = chrom if chrom.startswith('chr') else f"chr{chrom}"
            ncbi_chrom = ucsc_to_ncbi.get(ucsc_chrom, ucsc_chrom)
            
            coord_key = f"{ncbi_chrom}:{pos}"
            return coord_to_rsid.get(coord_key, coord_key) # Fallback to coordinate identifier if variant isn't inside dbSNP
        return clean_id

    snp_final.index = [convert_to_rsid(idx) for idx in snp_final.index]
    snp_final = snp_final[~snp_final.index.duplicated(keep='first')]

    # 6. COVARIATE MERGE
    cov_full = pd.read_csv("$COV_FILE", sep='\t')
    cov = cov_full[['externalsampleid'] + [c for c in cov_full.columns if c not in meta.columns]].drop_duplicates(subset='externalsampleid')
    meta_aligned = meta_unique[meta_unique['externalsubjectid'].isin(final_aligned_subjects)].copy()
    cov_merged = pd.merge(pd.merge(meta_aligned, cov, on='externalsampleid'), pca, left_on='externalsubjectid', right_on='IID')
    cov_merged = cov_merged.sort_values('RIN_score', ascending=False).drop_duplicates(subset='externalsubjectid')
    cov_merged['is_als'] = cov_merged['subject_group'].apply(lambda x: 1 if 'ALS' in str(x) else 0)
    cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map({'male': 1, 'female': 0})
    cov_final = cov_merged.set_index('externalsubjectid').reindex(final_aligned_subjects)[['sex_bin', 'age_at_death', 'is_als', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']].T

    # 7. SAVE OUTPUTS
    splicing_final.to_csv("$OUTDIR/splicing_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="geneid")
    snp_final.to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt", sep='\t', index=True, index_label="snpid")
    cov_final.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t', index=True, index_label="id", quoting=0)
    print(f"SUCCESS: Generated fully filtered datasets for {len(final_aligned_subjects)} samples.")
except Exception as e:
    import traceback; traceback.print_exc(); sys.exit(1)
EOF

# ─────────────────────────────────────────────────────────────────
# STEP 7: Coordinate Mapping Location Files
# ─────────────────────────────────────────────────────────────────
echo "Generating deduplicated location files for $TISSUE..."

# Generate SNP Location file tracking back to the dynamic rsIDs while outputting standard RefSeq NC_ structure
echo -e "snpid\tchr\tpos" > "$OUTDIR/snp_location.txt"
awk -v ncbi_map="$MAP_FILE" -v rsid_map="$TMP_RSID_MAP" '
BEGIN {
    OFS="\t"
    while ((getline < ncbi_map) > 0) { n_map[$2] = $1 }
    close(ncbi_map)
    while ((getline < rsid_map) > 0) { r_map[$1] = $2 }
    close(rsid_map)
}
{
    ucsc = ($1 ~ /^chr/) ? $1 : "chr"$1;
    ncbi = (ucsc in n_map) ? n_map[ucsc] : ucsc;
    
    coord_key = ncbi":"$4;
    final_id = (coord_key in r_map) ? r_map[coord_key] : coord_key;
    
    print final_id, ncbi, $4;
}' "$BIM" | sort -u -k1,1 >> "$OUTDIR/snp_location.txt"

# Splicing target location map (Tracks structural leafcutter BED data)
echo -e "geneid\tchr\tleft\tright" > "$OUTDIR/splicing_location.txt"
awk 'BEGIN{OFS="\t"} { print $4, $1, $2, $3 }' "$SPLICING_LOC_SRC" | sort -u -k1,1 >> "$OUTDIR/splicing_location.txt"

rm -f "$SPLICING_LOC_SRC" "$TMP_RSID_MAP"
echo "Location files complete for $TISSUE."

# ─────────────────────────────────────────────────────────────────
# STEP 8: Final Summary Output
# ─────────────────────────────────────────────────────────────────
echo "============================================"
echo "  sQTL Prep complete for $TISSUE"
echo "  Output directory: $OUTDIR"
echo "  $(date)"
echo "============================================"

echo "Output files:"
ls -lh "$OUTDIR"/*.txt 2>/dev/null || echo "No output files found"
