#!/bin/bash
#SBATCH --job-name=run_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/run_eQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/run_eQTL_%j.err
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
module load R

# ── Tissue argument ───────────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "ERROR: No tissue specified."
    echo "Usage: sbatch run_eQTL.sh \"Cerebellum\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  Matrix eQTL run for tissue : $TISSUE"
echo "  SLURM job ID               : ${SLURM_JOB_ID:-local}"
echo "  $(date)"
echo "============================================"

# ── Paths ─────────────────────────────────────────────────────────────
BASE=/home/zw529/donglab/data/target_ALS
PIPELINE=/home/zw529/donglab/pipelines/scripts/eQTL

INDIR=$BASE/$TISSUE_DIR/eQTL
OUTDIR=$INDIR/results
mkdir -p $OUTDIR

SNP_FILE=$INDIR/snp_${TISSUE_DIR}.txt
EXPR_FILE=$INDIR/expression_${TISSUE_DIR}.txt
COV_FILE=$INDIR/covariates_${TISSUE_DIR}_encoded.txt
GENE_LOC=$INDIR/gene_location.txt
SNP_LOC=$INDIR/snp_location.txt
OUTPUT_PREFIX=$OUTDIR/${TISSUE_DIR}_eQTL

# ── Pre-flight: check input files exist ───────────────────────────────
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
    echo "Run: sbatch prep_eQTL.sh \"$TISSUE\" first."
    exit 1
fi

# ── Pre-flight: verify sample counts match ────────────────────────────
echo ""
echo "[0b] Verifying sample count alignment..."

python3 << EOF
import sys

def col_count(fpath):
    with open(fpath) as f:
        return len(f.readline().rstrip('\n').split('\t')) - 1

n_snp  = col_count("$SNP_FILE")
n_expr = col_count("$EXPR_FILE")
n_cov  = col_count("$COV_FILE")

print(f"  SNP samples        : {n_snp}",  flush=True)
print(f"  Expression samples : {n_expr}", flush=True)
print(f"  Covariate samples  : {n_cov}",  flush=True)

ok = n_snp == n_expr == n_cov
if ok:
    print(f"  All sample counts match: {n_snp}", flush=True)
else:
    if n_expr != n_snp:
        print(f"  ERROR: Expression ({n_expr}) != SNP ({n_snp})", flush=True)
    if n_cov != n_snp:
        print(f"  ERROR: Covariates ({n_cov}) != SNP ({n_snp})", flush=True)
    print(f"  Re-run prep_eQTL.sh to fix alignment.", flush=True)
    sys.exit(1)
EOF

# ── Run Matrix eQTL ───────────────────────────────────────────────────
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

# ── Post-run summary ──────────────────────────────────────────────────
echo ""
echo "[2] Results summary for $TISSUE..."

CIS_FILE="${OUTPUT_PREFIX}.cis.txt"
TRANS_FILE="${OUTPUT_PREFIX}.trans.txt"
PDF_FILE="${OUTPUT_PREFIX}.pdf"

if [ -f "$CIS_FILE" ]; then
    N_CIS=$(tail -n +2 "$CIS_FILE" | wc -l)
    echo "  Cis eQTLs detected  : $N_CIS"
    echo "  Top 5 cis eQTLs:"
    head -6 "$CIS_FILE" | column -t
else
    echo "  WARNING: Cis output not found: $CIS_FILE"
fi

if [ -f "$TRANS_FILE" ]; then
    N_TRANS=$(tail -n +2 "$TRANS_FILE" | wc -l)
    echo "  Trans eQTLs detected: $N_TRANS"
    echo "  Top 5 trans eQTLs:"
    head -6 "$TRANS_FILE" | column -t
else
    echo "  WARNING: Trans output not found: $TRANS_FILE"
fi

[ -f "$PDF_FILE" ] && echo "  QQ plot saved to    : $PDF_FILE"

echo ""
echo "============================================"
echo "  Matrix eQTL complete for : $TISSUE"
echo "  Results in               : $OUTDIR"
echo "  $(date)"
echo "============================================"
