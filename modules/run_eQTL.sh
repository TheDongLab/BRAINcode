#!/bin/bash
#SBATCH --job-name=run_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/run_eQTL_%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/run_eQTL_%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

###########################################
# run_eQTL.sh
# Purpose: Run Matrix eQTL (cis + trans) for a given tissue type.
#           Expects prep_eQTL.sh to have been run first.
#
# Usage:
#   sbatch run_eQTL.sh "Cerebellum"
#   sbatch run_eQTL.sh "Frontal Cortex"
#   bash   run_eQTL.sh "Cerebellum"
###########################################

set -euo pipefail

# ---- Tissue argument ------------------------------------------------
if [ $# -lt 1 ]; then
    echo "ERROR: No tissue specified."
    echo "Usage: sbatch run_eQTL.sh \"Cerebellum\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  Matrix eQTL run for tissue: $TISSUE"
echo "  $(date)"
echo "============================================"

# ---- Paths ----------------------------------------------------------
BASE=/home/zw529/donglab/data/target_ALS
PIPELINE=/home/zw529/donglab/pipelines/scripts/eQTL

INDIR=$BASE/$TISSUE_DIR/eQTL
OUTDIR=$INDIR/results
mkdir -p $OUTDIR

SNP_FILE=$INDIR/snp_${TISSUE_DIR}.txt
EXPR_FILE=$INDIR/expression_${TISSUE_DIR}.txt
COV_FILE=$INDIR/covariates_${TISSUE_DIR}.txt
GENE_LOC=$INDIR/gene_location.txt
SNP_LOC=$INDIR/snp_location.txt
OUTPUT_PREFIX=$OUTDIR/${TISSUE_DIR}_eQTL

# ---- Pre-flight checks ----------------------------------------------
echo ""
echo "[0] Checking input files exist..."

MISSING=0
for f in "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" "$GENE_LOC" "$SNP_LOC"; do
    if [ ! -f "$f" ]; then
        echo "  MISSING: $f"
        MISSING=1
    else
        SIZE=$(du -sh "$f" | cut -f1)
        echo "  OK ($SIZE): $f"
    fi
done

if [ "$MISSING" -eq 1 ]; then
    echo ""
    echo "ERROR: One or more input files are missing."
    echo "Run prep_eQTL.sh \"$TISSUE\" first."
    exit 1
fi

# ---- Verify sample ID alignment across files ------------------------
echo ""
echo "[0b] Verifying sample ID alignment across expression, SNP, and covariate files..."

python3 << EOF
import sys

expr_file = "$EXPR_FILE"
snp_file  = "$SNP_FILE"
cov_file  = "$COV_FILE"

def get_header_ids(fpath, delimiter='\t', skip_first=True):
    with open(fpath) as f:
        header = f.readline().rstrip('\n').split(delimiter)
    return header[1:] if skip_first else header   # skip row-label column

expr_ids = get_header_ids(expr_file, '\t')
snp_ids  = get_header_ids(snp_file, '\t')

# Covariates: detect orientation
with open(cov_file) as f:
    cov_header = f.readline().rstrip('\n').split('\t')
cov_ids = cov_header[1:]   # assume wide format; adjust if long

expr_set = set(expr_ids)
snp_set  = set(snp_ids)
cov_set  = set(cov_ids)

print(f"  Expression samples : {len(expr_ids)}", flush=True)
print(f"  SNP samples        : {len(snp_ids)}", flush=True)
print(f"  Covariate samples  : {len(cov_ids)}", flush=True)

# Cross-check: expression vs SNP (these must match — different ID types)
# They use different ID systems (HRA vs HDA) so direct overlap won't be full
# but the ORDER must be consistent — we rely on prep_eQTL.sh having aligned them
# via the subject-level join. Warn if counts differ significantly.
if len(expr_ids) != len(snp_ids):
    print(f"  WARNING: Expression has {len(expr_ids)} samples but SNP matrix has {len(snp_ids)} samples.", flush=True)
    print(f"  Matrix eQTL requires matched columns — check prep_eQTL.sh output.", flush=True)
    sys.exit(1)
else:
    print(f"  Sample counts match: {len(expr_ids)} samples in both expression and SNP files.", flush=True)

if len(cov_ids) > 0 and len(cov_ids) != len(expr_ids):
    print(f"  WARNING: Covariate file has {len(cov_ids)} samples vs {len(expr_ids)} in expression.", flush=True)
else:
    print(f"  Covariate sample count OK.", flush=True)
EOF

# ---- Run Matrix eQTL ------------------------------------------------
echo ""
echo "[1] Running Matrix eQTL (cis + trans)..."
echo "  SNP file    : $SNP_FILE"
echo "  Expr file   : $EXPR_FILE"
echo "  Cov file    : $COV_FILE"
echo "  Gene loc    : $GENE_LOC"
echo "  SNP loc     : $SNP_LOC"
echo "  Output      : $OUTPUT_PREFIX"
echo ""

Rscript $PIPELINE/_eQTL.R \
    "$SNP_FILE" \
    "$EXPR_FILE" \
    "$COV_FILE" \
    "$OUTPUT_PREFIX" \
    "$GENE_LOC" \
    "$SNP_LOC"

# ---- Post-run summary -----------------------------------------------
echo ""
echo "[2] Results summary for $TISSUE..."

CIS_FILE="${OUTPUT_PREFIX}.cis.txt"
TRANS_FILE="${OUTPUT_PREFIX}.trans.txt"
PDF_FILE="${OUTPUT_PREFIX}.pdf"

if [ -f "$CIS_FILE" ]; then
    N_CIS=$(tail -n +2 $CIS_FILE | wc -l)
    echo "  Cis eQTLs detected  : $N_CIS"
    echo "  Top 5 cis eQTLs:"
    head -6 $CIS_FILE | column -t
else
    echo "  WARNING: Cis output file not found: $CIS_FILE"
fi

if [ -f "$TRANS_FILE" ]; then
    N_TRANS=$(tail -n +2 $TRANS_FILE | wc -l)
    echo "  Trans eQTLs detected: $N_TRANS"
    echo "  Top 5 trans eQTLs:"
    head -6 $TRANS_FILE | column -t
else
    echo "  WARNING: Trans output file not found: $TRANS_FILE"
fi

if [ -f "$PDF_FILE" ]; then
    echo "  QQ plot saved to    : $PDF_FILE"
fi

echo ""
echo "============================================"
echo "  Matrix eQTL complete for: $TISSUE"
echo "  Results in: $OUTDIR"
echo "  $(date)"
echo "============================================"
