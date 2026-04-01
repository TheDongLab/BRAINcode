#!/bin/bash
#SBATCH --job-name=prep_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=499G

###########################################
# prep_eQTL.sh
# Purpose: Subset all Matrix eQTL inputs to a given tissue type,
#           and prepare gene/SNP location files.
#
# Logic: Keeps subjects who have BOTH:
#   - RNAseq for the specified tissue in expression_sample_metadata.csv
#   - WGS data (any tissue) in wgs_samples_for_vcf_merge.csv
#
# Usage:
#   sbatch prep_eQTL.sh "Cerebellum"
#   sbatch prep_eQTL.sh "Frontal Cortex"
#   bash   prep_eQTL.sh "Cerebellum"
#
# Output: $BASE/[tissue_dir]/eQTL/
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
# Convert spaces to underscores for use in directory/file names
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  eQTL Prep for tissue: $TISSUE"
echo "  Directory label     : $TISSUE_DIR"
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

# Temp files namespaced by tissue to allow parallel runs
TMP_META=/tmp/${TISSUE_DIR}_meta.txt
TMP_WGS=/tmp/${TISSUE_DIR}_wgs_map.txt
TMP_MATCHED=/tmp/${TISSUE_DIR}_matched.txt
TMP_HRA=/tmp/${TISSUE_DIR}_HRA_ids.txt
TMP_HRA_UNDER=/tmp/${TISSUE_DIR}_HRA_ids_underscore.txt
TMP_HDA=/tmp/${TISSUE_DIR}_HDA_ids.txt

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
    echo "Check tissue label matches exactly — available labels:"
    tail -n +2 $META | awk -F',' '{print $3}' | sort | uniq -c | sort -rn
    exit 1
fi

tail -n +2 $WGS_MAP \
    | awk -F',' '{print $1","$2}' \
    | sort -t',' -k1,1 \
    > $TMP_WGS

N_WGS=$(wc -l < $TMP_WGS)
echo "  Subjects with WGS data             : $N_WGS"

join -t',' -1 1 -2 1 \
    $TMP_META \
    $TMP_WGS \
    > $TMP_MATCHED

N_MATCHED=$(wc -l < $TMP_MATCHED)
echo "  Subjects with $TISSUE RNA + WGS    : $N_MATCHED"

if [ "$N_MATCHED" -eq 0 ]; then
    echo "ERROR: No matched samples found. Check subjectid formatting between metadata and WGS map."
    exit 1
fi

awk -F',' '{print $2}' $TMP_MATCHED > $TMP_HRA
awk -F',' '{print $3}' $TMP_MATCHED > $TMP_HDA

echo "  Sample HRA IDs (expression): $(head -3 $TMP_HRA | tr '\n' ' ')..."
echo "  Sample HDA IDs (SNP):        $(head -3 $TMP_HDA | tr '\n' ' ')..."

##############################################
# STEP 2: Subset expression matrix
##############################################
echo ""
echo "[2] Subsetting expression matrix..."

sed 's/-/_/g' $TMP_HRA > $TMP_HRA_UNDER

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
# STEP 3: Subset covariates
##############################################
echo ""
echo "[3] Subsetting covariates..."

python3 << EOF
ids_file = "$TMP_HRA"
cov_file = "$COV"
out_file = "$OUTDIR/covariates_${TISSUE_DIR}.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

with open(cov_file) as f:
    header = f.readline().rstrip('\n').split('\t')

sample_in_header = any(h in keep_ids for h in header)

with open(cov_file) as f, open(out_file, 'w') as out:
    header = f.readline().rstrip('\n').split('\t')
    if sample_in_header:
        keep_idx = [0] + [i for i, h in enumerate(header) if h in keep_ids]
        print(f"  Covariate format  : WIDE (samples as columns)", flush=True)
        print(f"  Keeping {len(keep_idx)-1} sample columns", flush=True)
        out.write('\t'.join(header[i] for i in keep_idx) + '\n')
        for line in f:
            fields = line.rstrip('\n').split('\t')
            out.write('\t'.join(fields[i] for i in keep_idx) + '\n')
    else:
        print(f"  Covariate format  : LONG (samples as rows)", flush=True)
        out.write('\t'.join(header) + '\n')
        n = 0
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if fields[0] in keep_ids:
                out.write(line)
                n += 1
        print(f"  Kept {n} sample rows", flush=True)

print(f"  Written to: {out_file}", flush=True)
EOF

echo "  Covariates done."

##############################################
# STEP 4: Subset and transpose SNP matrix
##############################################
echo ""
echo "[4] Subsetting and transposing SNP matrix (chunked)..."

python3 << EOF
import re

ids_file = "$TMP_HDA"
raw_file = "$RAW"
out_file = "$OUTDIR/snp_${TISSUE_DIR}.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

print(f"  HDA IDs to keep: {len(keep_ids)}", flush=True)

def normalize_hda(iid):
    return re.sub(r'-b\d+\$', '', iid)

rows = []
with open(raw_file) as f:
    header = f.readline().split()
    snp_names = header[6:]
    for line in f:
        fields = line.split()
        iid_raw = fields[1]
        iid_norm = normalize_hda(iid_raw)
        matched_id = iid_norm if iid_norm in keep_ids else (iid_raw if iid_raw in keep_ids else None)
        if matched_id:
            rows.append((matched_id, fields[6:]))

print(f"  Samples found in SNP matrix : {len(rows)}", flush=True)
missing = keep_ids - set(r[0] for r in rows)
if missing:
    print(f"  WARNING: {len(missing)} HDA IDs not found in SNP matrix", flush=True)
    print(f"  First 5 missing: {list(missing)[:5]}", flush=True)

sample_ids = [r[0] for r in rows]
n_samples   = len(rows)
n_snps      = len(snp_names)
CHUNK       = 50000

with open(out_file, 'w') as out:
    out.write('snpid\t' + '\t'.join(sample_ids) + '\n')
    for chunk_start in range(0, n_snps, CHUNK):
        chunk_end = min(chunk_start + CHUNK, n_snps)
        for j in range(chunk_start, chunk_end):
            snp_row = [rows[i][1][j] for i in range(n_samples)]
            out.write(snp_names[j] + '\t' + '\t'.join(snp_row) + '\n')
        print(f"  Written SNPs {chunk_start}–{chunk_end} of {n_snps}", flush=True)

print(f"  Done: {n_snps} SNPs x {n_samples} samples -> {out_file}", flush=True)
EOF

echo "  SNP matrix done."

##############################################
# STEP 5: SNP location file from .bim
##############################################
echo ""
echo "[5] Generating SNP location file..."

SNP_LOC=$OUTDIR/snp_location.txt
echo -e "snpid\tchr\tpos" > $SNP_LOC
awk 'BEGIN{OFS="\t"} {print $2, "chr"$1, $4}' $BIM >> $SNP_LOC

N_SNPS=$(tail -n +2 $SNP_LOC | wc -l)
echo "  SNPs in location file : $N_SNPS"
echo "  Written to: $SNP_LOC"

##############################################
# STEP 6: Gene location file from gencode v49 BED6
##############################################
echo ""
echo "[6] Generating gene location file..."

GENE_LOC=$OUTDIR/gene_location.txt
echo -e "geneid\tchr\tleft\tright" > $GENE_LOC
awk 'BEGIN{OFS="\t"} {
    gene=$4; sub(/\.[0-9]+$/, "", gene)
    print gene, $1, $2, $3
}' $GTF_BED6 >> $GENE_LOC

N_GENES=$(tail -n +2 $GENE_LOC | wc -l)
echo "  Genes in location file : $N_GENES"
echo "  Written to: $GENE_LOC"

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
print(f"  Expression genes              : {len(expr_genes)}", flush=True)
print(f"  Found in gene location file   : {found}", flush=True)
print(f"  Missing from location file    : {len(missing)}", flush=True)
if missing:
    print(f"  First 5 missing gene IDs      : {missing[:5]}", flush=True)
EOF

##############################################
# STEP 7: Summary
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
