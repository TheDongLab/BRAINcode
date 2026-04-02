#!/bin/bash
#SBATCH --job-name=prep_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=499G

###########################################
# prep_eQTL.sh
# Purpose: Subset all Matrix eQTL inputs to a given tissue type
#           and prepare gene/SNP location files.
#
# Logic: Keeps subjects with BOTH:
#   - RNAseq for the specified tissue in expression_sample_metadata.csv
#   - WGS data (any tissue) in wgs_samples_for_vcf_merge.csv
#
# Covariate file is LONG format (samples as rows, covariates as columns)
# and is transposed to WIDE format for Matrix eQTL.
#
# A final alignment step enforces that expression, SNP, and covariate
# files contain exactly the same samples in the same order.
#
# Usage:
#   sbatch prep_eQTL.sh "Cerebellum"
#   sbatch prep_eQTL.sh "Frontal Cortex"
#   bash   prep_eQTL.sh "Cerebellum"
#
# Output: $BASE/[TISSUE_DIR]/eQTL/
###########################################

set -euo pipefail

# ---- Tissue argument ------------------------------------------------
if [ $# -lt 1 ]; then
    echo "ERROR: No tissue specified."
    echo "Usage: sbatch prep_eQTL.sh \"Cerebellum\""
    echo "       sbatch prep_eQTL.sh \"Frontal Cortex\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  eQTL Prep for tissue : $TISSUE"
echo "  Directory label      : $TISSUE_DIR"
echo "  $(date)"
echo "============================================"

# ---- Paths ----------------------------------------------------------
BASE=/home/zw529/donglab/data/target_ALS
PLINK=$BASE/eQTL/plink
REFS=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode

EXPR=$BASE/eQTL/expression_matrix.txt
RAW=$PLINK/joint_autosomes_matrixEQTL.raw
META=$BASE/eQTL/expression_sample_metadata.csv
WGS_MAP=$BASE/eQTL/wgs_samples_for_vcf_merge.csv
COV=$BASE/eQTL/covariates.tsv
BIM=$PLINK/joint_autosomes_filtered_bed.bim
GTF_BED6=$REFS/gencode.v49.annotation.gene.bed6

OUTDIR=$BASE/$TISSUE_DIR/eQTL
mkdir -p $OUTDIR

# Temp files written to OUTDIR (not /tmp) so they persist across compute nodes
TMP_META=$OUTDIR/tmp_meta.txt
TMP_WGS=$OUTDIR/tmp_wgs_map.txt
TMP_MATCHED=$OUTDIR/tmp_matched.txt
TMP_HRA=$OUTDIR/tmp_HRA_ids.txt
TMP_HRA_UNDER=$OUTDIR/tmp_HRA_ids_underscore.txt
TMP_HDA=$OUTDIR/tmp_HDA_ids.txt
TMP_FINAL_SAMPLES=$OUTDIR/tmp_final_sample_list.txt   # ordered intersection after alignment

##############################################
# STEP 1: Build sample mapping table
##############################################
echo ""
echo "[1] Building sample mapping (HRA <-> subjectid <-> HDA) for: $TISSUE"

tail -n +2 $META \
    | awk -F',' -v tissue="$TISSUE" '$3 == tissue {print $2","$1}' \
    | sort -t',' -k1,1 \
    > $TMP_META

N_META=$(wc -l < $TMP_META)
echo "  RNA samples in metadata for $TISSUE : $N_META"

if [ "$N_META" -eq 0 ]; then
    echo "ERROR: No samples found for tissue '$TISSUE'."
    echo "Check tissue label exactly — available labels:"
    tail -n +2 $META | awk -F',' '{print $3}' | sort | uniq -c | sort -rn
    exit 1
fi

tail -n +2 $WGS_MAP \
    | awk -F',' '{print $1","$2}' \
    | sort -t',' -k1,1 \
    > $TMP_WGS

N_WGS=$(wc -l < $TMP_WGS)
echo "  Subjects with WGS data             : $N_WGS"

# Inner join on subjectid → subjectid, HRA_id, HDA_id
join -t',' -1 1 -2 1 \
    $TMP_META \
    $TMP_WGS \
    > $TMP_MATCHED

N_MATCHED=$(wc -l < $TMP_MATCHED)
echo "  Subjects with $TISSUE RNA + WGS    : $N_MATCHED"

if [ "$N_MATCHED" -eq 0 ]; then
    echo "ERROR: No matched samples found."
    exit 1
fi

awk -F',' '{print $2}' $TMP_MATCHED > $TMP_HRA
awk -F',' '{print $3}' $TMP_MATCHED > $TMP_HDA
sed 's/-/_/g' $TMP_HRA > $TMP_HRA_UNDER

echo "  Sample HRA IDs (expression): $(head -3 $TMP_HRA | tr '\n' ' ')..."
echo "  Sample HDA IDs (SNP):        $(head -3 $TMP_HDA | tr '\n' ' ')..."

##############################################
# STEP 2: Subset expression matrix
# Expression matrix columns use underscores (CGND_HRA_03594)
##############################################
echo ""
echo "[2] Subsetting expression matrix..."

python3 << EOF
ids_file  = "$TMP_HRA_UNDER"
expr_file = "$EXPR"
out_file  = "$OUTDIR/expression_${TISSUE_DIR}.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

with open(expr_file) as f, open(out_file, 'w') as out:
    header = f.readline().rstrip('\n').split('\t')
    keep_idx = [0] + [i for i, h in enumerate(header) if h in keep_ids]
    found_ids = [header[i] for i in keep_idx[1:]]
    missing = keep_ids - set(found_ids)

    print(f"  Requested samples : {len(keep_ids)}", flush=True)
    print(f"  Found in matrix   : {len(found_ids)}", flush=True)
    if missing:
        print(f"  WARNING - {len(missing)} IDs not found in expression matrix", flush=True)
        print(f"  First 5 missing   : {list(missing)[:5]}", flush=True)

    out.write('\t'.join(header[i] for i in keep_idx) + '\n')
    for line in f:
        fields = line.rstrip('\n').split('\t')
        out.write('\t'.join(fields[i] for i in keep_idx) + '\n')

print(f"  Written to: {out_file}", flush=True)
EOF

echo "  Expression matrix done."

##############################################
# STEP 3: Subset and transpose SNP matrix
# .raw samples are ROWS; Matrix eQTL wants SNPs as ROWS
# Strip -b38 suffix from IIDs before matching
##############################################
echo ""
echo "[3] Subsetting and transposing SNP matrix (chunked)..."

python3 << EOF
import re

ids_file = "$TMP_HDA"
raw_file = "$RAW"
out_file = "$OUTDIR/snp_${TISSUE_DIR}.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

print(f"  HDA IDs to keep: {len(keep_ids)}", flush=True)

def normalize_hda(iid):
    return re.sub(r'-b\d+$', '', iid)

rows = []
with open(raw_file) as f:
    header = f.readline().split()
    snp_names = header[6:]
    for line in f:
        fields = line.split()
        iid_raw   = fields[1]
        iid_norm  = normalize_hda(iid_raw)
        matched   = iid_norm if iid_norm in keep_ids else (iid_raw if iid_raw in keep_ids else None)
        if matched:
            rows.append((matched, fields[6:]))

print(f"  Samples found in SNP matrix : {len(rows)}", flush=True)
missing = keep_ids - set(r[0] for r in rows)
if missing:
    print(f"  WARNING: {len(missing)} HDA IDs not found in SNP matrix", flush=True)
    print(f"  First 5 missing: {list(missing)[:5]}", flush=True)

sample_ids = [r[0] for r in rows]
n_samples  = len(rows)
n_snps     = len(snp_names)
CHUNK      = 50000

with open(out_file, 'w') as out:
    out.write('snpid\t' + '\t'.join(sample_ids) + '\n')
    for chunk_start in range(0, n_snps, CHUNK):
        chunk_end = min(chunk_start + CHUNK, n_snps)
        for j in range(chunk_start, chunk_end):
            snp_row = [rows[i][1][j] for i in range(n_samples)]
            out.write(snp_names[j] + '\t' + '\t'.join(snp_row) + '\n')
        print(f"  Written SNPs {chunk_start}-{chunk_end} of {n_snps}", flush=True)

print(f"  Done: {n_snps} SNPs x {n_samples} samples -> {out_file}", flush=True)
EOF

echo "  SNP matrix done."

##############################################
# STEP 4: Subset covariates
# LONG format: rows = samples, cols = covariate variables
# We subset rows to tissue samples, then TRANSPOSE to wide format for Matrix eQTL (rows = covariates, cols = samples)
# Numeric covariates only — sex/ethnicity etc need to be encoded
##############################################
echo ""
echo "[4] Subsetting and transposing covariates..."

python3 << EOF
ids_file = "$TMP_HRA"    # hyphens, matching covariates.tsv
cov_file = "$COV"
out_file = "$OUTDIR/covariates_${TISSUE_DIR}.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

# Read long-format covariate file (rows = samples)
header = None
kept_rows = {}   # sample_id -> {covariate: value}

with open(cov_file) as f:
    header = f.readline().rstrip('\n').split('\t')
    sample_col = header[0]   # 'externalsampleid'
    cov_cols   = header[1:]  # covariate column names

    for line in f:
        fields = line.rstrip('\n').split('\t')
        sid = fields[0]
        if sid in keep_ids:
            kept_rows[sid] = fields[1:]

found = len(kept_rows)
missing = keep_ids - set(kept_rows.keys())
print(f"  Requested samples : {len(keep_ids)}", flush=True)
print(f"  Found in cov file : {found}", flush=True)
if missing:
    print(f"  WARNING - {len(missing)} IDs not in covariates", flush=True)
    print(f"  First 5 missing   : {list(missing)[:5]}", flush=True)

# Transpose: write wide format — rows = covariates, cols = samples
# Header row: covariate_name  sample1  sample2 ...
sample_order = sorted(kept_rows.keys())   # consistent ordering

with open(out_file, 'w') as out:
    out.write('covariate\t' + '\t'.join(sample_order) + '\n')
    for i, cov_name in enumerate(cov_cols):
        row_vals = [kept_rows[s][i] if s in kept_rows else 'NA' for s in sample_order]
        out.write(cov_name + '\t' + '\t'.join(row_vals) + '\n')

print(f"  Transposed {len(cov_cols)} covariates x {len(sample_order)} samples", flush=True)
print(f"  Written to: {out_file}", flush=True)
EOF

echo "  Covariates done."

##############################################
# STEP 5: Align all three files to exact same sample set and order
#
# The expression matrix uses HRA IDs (underscores)
# The SNP matrix uses HDA IDs
# The covariate matrix uses HRA IDs (hyphens)
#
# Matrix eQTL matches samples by COLUMN POSITION not by name across files, so expression col 1 must correspond to SNP col 1 etc.
# We use the matched table (subjectid, HRA, HDA) to enforce this.
##############################################
echo ""
echo "[5] Aligning samples across expression, SNP, and covariate files..."

python3 << EOF
import re

matched_file = "$TMP_MATCHED"      # subjectid, HRA_id, HDA_id
expr_file    = "$OUTDIR/expression_${TISSUE_DIR}.txt"
snp_file     = "$OUTDIR/snp_${TISSUE_DIR}.txt"
cov_file     = "$OUTDIR/covariates_${TISSUE_DIR}.txt"
final_list   = "$TMP_FINAL_SAMPLES"   # written for logging

# Load the matched table: subjectid -> (HRA_underscore, HDA)
subject_to_hra = {}
subject_to_hda = {}
with open(matched_file) as f:
    for line in f:
        parts = line.strip().split(',')
        subj, hra, hda = parts[0], parts[1], parts[2]
        hra_under = hra.replace('-', '_')
        subject_to_hra[subj] = hra_under
        subject_to_hda[subj] = hda

def normalize_hda(iid):
    return re.sub(r'-b\d+$', '', iid)

# Get samples actually present in each file
def get_col_ids(fpath):
    with open(fpath) as f:
        return f.readline().rstrip('\n').split('\t')[1:]

expr_ids_present = set(get_col_ids(expr_file))
snp_ids_raw      = get_col_ids(snp_file)
snp_ids_present  = set(normalize_hda(x) for x in snp_ids_raw)

# For covariates: column headers are sample IDs (HRA with hyphens)
cov_ids_raw     = get_col_ids(cov_file)
cov_ids_present = set(cov_ids_raw)

# Find subjects present in ALL three files
final_subjects = []
for subj, hra_under in subject_to_hra.items():
    hda = normalize_hda(subject_to_hda[subj])
    hra_hyphen = hra_under.replace('_', '-')
    if hra_under in expr_ids_present and hda in snp_ids_present and hra_hyphen in cov_ids_present:
        final_subjects.append(subj)

print(f"  Subjects in matched table      : {len(subject_to_hra)}", flush=True)
print(f"  Subjects present in expr       : {sum(1 for s,h in subject_to_hra.items() if h in expr_ids_present)}", flush=True)
print(f"  Subjects present in SNP        : {sum(1 for s in subject_to_hra if normalize_hda(subject_to_hda[s]) in snp_ids_present)}", flush=True)
print(f"  Subjects present in covariates : {sum(1 for s,h in subject_to_hra.items() if h.replace('_','-') in cov_ids_present)}", flush=True)
print(f"  Final aligned sample count     : {len(final_subjects)}", flush=True)

# Write ordered HRA (underscore) and HDA lists for re-subsetting
hra_ordered = [subject_to_hra[s] for s in final_subjects]
hda_ordered = [normalize_hda(subject_to_hda[s]) for s in final_subjects]
hra_hyphen_ordered = [h.replace('_', '-') for h in hra_ordered]

with open(final_list, 'w') as f:
    for subj, hra, hda in zip(final_subjects, hra_ordered, hda_ordered):
        f.write(f"{subj}\t{hra}\t{hda}\n")

# --- Re-subset expression to final ordered samples ---
print(f"  Re-subsetting expression matrix...", flush=True)
with open(expr_file) as f:
    header = f.readline().rstrip('\n').split('\t')
    idx_map = {h: i for i, h in enumerate(header)}
    keep_idx = [0] + [idx_map[h] for h in hra_ordered if h in idx_map]
    rows = ['\t'.join(header[i] for i in keep_idx)]
    for line in f:
        fields = line.rstrip('\n').split('\t')
        rows.append('\t'.join(fields[i] for i in keep_idx))
with open(expr_file, 'w') as f:
    f.write('\n'.join(rows) + '\n')
print(f"  Expression re-subset: {len(keep_idx)-1} samples", flush=True)

# --- Re-subset SNP matrix to final ordered samples ---
# SNP matrix is huge — stream it
print(f"  Re-subsetting SNP matrix (streaming)...", flush=True)
snp_tmp = snp_file + ".tmp"
hda_set_ordered = hda_ordered   # list preserving order

with open(snp_file) as fin, open(snp_tmp, 'w') as fout:
    header = fin.readline().rstrip('\n').split('\t')
    # Build index map from normalized HDA -> column index
    norm_to_idx = {}
    for i, h in enumerate(header[1:], start=1):
        norm_to_idx[normalize_hda(h)] = i
    keep_idx = [0] + [norm_to_idx[h] for h in hda_ordered if h in norm_to_idx]
    fout.write('\t'.join(header[i] for i in keep_idx) + '\n')
    for line in fin:
        fields = line.rstrip('\n').split('\t')
        fout.write('\t'.join(fields[i] for i in keep_idx) + '\n')

import os
os.replace(snp_tmp, snp_file)
print(f"  SNP re-subset: {len(keep_idx)-1} samples", flush=True)

# --- Re-subset covariate matrix to final ordered samples ---
print(f"  Re-subsetting covariate matrix...", flush=True)
with open(cov_file) as f:
    header = f.readline().rstrip('\n').split('\t')
    idx_map = {h: i for i, h in enumerate(header)}
    keep_idx = [0] + [idx_map[h] for h in hra_hyphen_ordered if h in idx_map]
    rows = ['\t'.join(header[i] for i in keep_idx)]
    for line in f:
        fields = line.rstrip('\n').split('\t')
        rows.append('\t'.join(fields[i] for i in keep_idx))
with open(cov_file, 'w') as f:
    f.write('\n'.join(rows) + '\n')
print(f"  Covariates re-subset: {len(keep_idx)-1} samples", flush=True)

print(f"  Alignment complete.", flush=True)
EOF

echo "  Alignment done."

##############################################
# STEP 6: SNP location file from .bim
##############################################
echo ""
echo "[6] Generating SNP location file..."

SNP_LOC=$OUTDIR/snp_location.txt
echo -e "snpid\tchr\tpos" > $SNP_LOC
awk 'BEGIN{OFS="\t"} {print $2, "chr"$1, $4}' $BIM >> $SNP_LOC

N_SNPS=$(tail -n +2 $SNP_LOC | wc -l)
echo "  SNPs in location file : $N_SNPS"

##############################################
# STEP 7: Gene location file from gencode v49 BED6
##############################################
echo ""
echo "[7] Generating gene location file..."

GENE_LOC=$OUTDIR/gene_location.txt
echo -e "geneid\tchr\tleft\tright" > $GENE_LOC

# extract only ENSG ID (before first underscore)
awk 'BEGIN{OFS="\t"} {
    gene=$4
    split(gene, a, "___")   # split on triple underscore
    gene_id=a[1]            # ENSGxxxxxx.version
    print gene_id, $1, $2, $3
}' $GTF_BED6 >> $GENE_LOC

N_GENES=$(tail -n +2 $GENE_LOC | wc -l)
echo "  Genes in location file : $N_GENES"

N_GENES=$(tail -n +2 $GENE_LOC | wc -l)
echo "  Genes in location file : $N_GENES"

python3 << EOF
expr_file = "$OUTDIR/expression_${TISSUE_DIR}.txt"
loc_file  = "$OUTDIR/gene_location.txt"

with open(loc_file) as f:
    f.readline()
    loc_genes = set(line.split('\t')[0] for line in f)
with open(expr_file) as f:
    f.readline()
    expr_genes = [line.split('\t')[0] for line in f]

found   = sum(1 for g in expr_genes if g in loc_genes)
missing = [g for g in expr_genes if g not in loc_genes]
print(f"  Expression genes            : {len(expr_genes)}", flush=True)
print(f"  Found in gene location file : {found}", flush=True)
print(f"  Missing                     : {len(missing)}", flush=True)
if missing:
    print(f"  First 5 missing             : {missing[:5]}", flush=True)
EOF

##############################################
# STEP 8: Cleanup temp files
##############################################
echo ""
echo "[8] Cleaning up temp files..."
rm -f $TMP_META $TMP_WGS $TMP_MATCHED $TMP_HRA $TMP_HRA_UNDER $TMP_HDA
echo "  Done. Final sample list retained at: $TMP_FINAL_SAMPLES"

##############################################
# STEP 9: Final summary
##############################################
echo ""
echo "============================================"
echo "  Prep complete for: $TISSUE"
echo "  Output directory : $OUTDIR"
echo ""
echo "  snp_${TISSUE_DIR}.txt         -> args[1] SNP_file_name"
echo "  expression_${TISSUE_DIR}.txt  -> args[2] expression_file_name"
echo "  covariates_${TISSUE_DIR}.txt  -> args[3] covariates_file_name"
echo "  gene_location.txt             -> args[5] gene_location_file_name"
echo "  snp_location.txt              -> args[6] snp_location_file_name"
echo ""
echo "  Next step: sbatch run_eQTL.sh \"$TISSUE\""
echo "  $(date)"
echo "============================================"
