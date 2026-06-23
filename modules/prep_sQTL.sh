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
                        out.write(f"{chr_part}\t{start}\t{end}\t{event_id}\n")
            print(f"STEP 0: Complete using encoding: {enc}", file=sys.stderr)
            break
        except UnicodeDecodeError: continue
except Exception as e:
    print(f"ERROR in STEP 0: {e}", file=sys.stderr); sys.exit(1)
EOF

##############################################
# STEP 1-6: Load, Filter, and Align Matrix Data
##############################################
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

    sub_to_hra, sub_to_raw_sample, final_aligned_subjects = {}, {}, []
    junc_counts = Counter()
    rep_dir_base = "/home/zw529/donglab/data/target_ALS/" + "$TISSUE_DIR" + "/RNAseq/Processed"

    for s in common_subjects:
        row = meta_unique[meta_unique['externalsubjectid'] == s]
        hra = row['hra_clean'].values[0]
        raw_sample = row['externalsampleid'].values[0] 
        matches = [c for c in splicing_headers if c.startswith(hra)]
        if matches:
            tsv_path = os.path.join(rep_dir_base, raw_sample, "leafcutter", "psi", f"{raw_sample}.leafcutter.PSI.tsv")
            if os.path.exists(tsv_path):
                final_aligned_subjects.append(s)
                sub_to_hra[s] = matches[0]
                try:
                    df_rep = pd.read_csv(tsv_path, sep=r'\s+', usecols=['chrom', 'start', 'end', 'reads'])
                    df_rep['j_id'] = df_rep['chrom'].astype(str) + ':' + df_rep['start'].astype(str) + '-' + df_rep['end'].astype(str)
                    junc_counts.update(df_rep.loc[df_rep['reads'] >= 5, 'j_id'].unique())
                except: pass

    if len(final_aligned_subjects) == 0:
        raise ValueError(f"CRITICAL ERROR: Zero sample files matched! Checked base: {rep_dir_base}")

    # Retain junctions with >=5 reads in >=10% of successfully mapped samples
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
    def clean_snp_id(fid):
        cid = fid.rsplit('_', 1)[0] if '_' in fid else fid
        return cid if str(cid).startswith('chr') else f"chr{cid.split(':')[0]}:{cid.split(':')[1]}" if ':' in str(cid) else cid
    snp_final.index = [clean_snp_id(idx) for idx in snp_final.index]
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

##############################################
# STEP 7: Coordinate Mapping Location Files
##############################################
echo "Generating deduplicated location files for $TISSUE..."
echo -e "snpid\tchr\tpos" > "$OUTDIR/snp_location.txt"
awk 'BEGIN{OFS="\t"} {chrom = ($1 ~ /^chr/) ? $1 : "chr"$1; print "chr"$1":"$4, chrom, $4;}' "$BIM" | sed 's/chrchr/chr/g' | sort -u -k1,1 >> "$OUTDIR/snp_location.txt"

echo -e "geneid\tchr\tleft\tright" > "$OUTDIR/splicing_location.txt"
awk 'BEGIN{OFS="\t"} { print $4, $1, $2, $3 }' "$SPLICING_LOC_SRC" | sort -u -k1,1 >> "$OUTDIR/splicing_location.txt"
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

# List output files for initial script diagnostics
echo "Output files:"
ls -lh "$OUTDIR"/*.txt 2>/dev/null || echo "No output files found"
