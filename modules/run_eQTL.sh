#!/bin/bash
#SBATCH --job-name=run_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/run_eQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/run_eQTL_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

###########################################
# run_eQTL.sh
# Purpose: Run Matrix eQTL (cis only) for a given tissue type,
#           then post-process results and generate figures.
#
# Expects prep_eQTL.sh to have been run first.
#
# Steps:
#   0.  Pre-flight checks (file existence + sample count alignment)
#   1.  Run Matrix eQTL via _eQTL.R
#   2.  Post-process: join coordinates, flag telomeric SNPs,
#       filter FDR, identify lead SNPs (_eQTL_postprocess.R)
#   3.  Manhattan plots by SNP position and gene position
#       (_eQTL_manhattan.R)
#   4.  Boxplots for top 100 gene-SNP pairs (_eQTL_boxplot.R)
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
PLINK=$BASE/eQTL/plink

INDIR=$BASE/$TISSUE_DIR/eQTL
OUTDIR=$INDIR/results
mkdir -p $OUTDIR

SNP_FILE=$INDIR/snp_${TISSUE_DIR}.txt
EXPR_FILE=$INDIR/expression_${TISSUE_DIR}.txt
COV_FILE=$INDIR/covariates_${TISSUE_DIR}_encoded.txt
GENE_LOC=$INDIR/gene_location.txt
SNP_LOC=$INDIR/snp_location.txt
BIM=$PLINK/joint_autosomes_filtered_bed.bim
OUTPUT_PREFIX=$OUTDIR/${TISSUE_DIR}_eQTL

CIS_FILE="${OUTPUT_PREFIX}.cis.txt"
FDR_THRESH=0.05
TOP_N=100

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

if n_snp == n_expr == n_cov:
    print(f"  All sample counts match: {n_snp}", flush=True)
else:
    if n_expr != n_snp:
        print(f"  ERROR: Expression ({n_expr}) != SNP ({n_snp})", flush=True)
    if n_cov != n_snp:
        print(f"  ERROR: Covariates ({n_cov}) != SNP ({n_snp})", flush=True)
    print(f"  Re-run prep_eQTL.sh to fix alignment.", flush=True)
    sys.exit(1)
EOF

# ── Step 1: Run Matrix eQTL ───────────────────────────────────────────
echo ""
echo "[1] Running Matrix eQTL (cis only)..."
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

# ── Step 2: Post-processing ───────────────────────────────────────────
echo ""
echo "[2] Post-processing results..."

if [ ! -f "$CIS_FILE" ]; then
    echo "  ERROR: Cis output not found: $CIS_FILE"
    echo "  Matrix eQTL may have failed — check the error log."
    exit 1
fi

N_CIS=$(tail -n +2 "$CIS_FILE" | wc -l)
echo "  Raw cis-eQTLs detected : $N_CIS"

Rscript $PIPELINE/_eQTL_postprocess.R \
    "$CIS_FILE" \
    "$SNP_LOC" \
    "$GENE_LOC" \
    "$OUTPUT_PREFIX" \
    "$FDR_THRESH" \
    "$TOP_N"

ANNOTATED_FILE="${OUTPUT_PREFIX}.full_annotated.txt"
LEAD_FILE="${OUTPUT_PREFIX}.lead_snps.txt"
FDR_FILE="${OUTPUT_PREFIX}.FDR${FDR_THRESH}.txt"
TOP_PAIRS="${OUTPUT_PREFIX}.top${TOP_N}_for_boxplot.txt"

N_FDR=$(tail -n +2 "$FDR_FILE" | wc -l)
N_LEAD=$(wc -l < "$LEAD_FILE")
echo "  FDR < $FDR_THRESH associations : $N_FDR"
echo "  Unique genes with lead SNP      : $N_LEAD"

# ── Step 3: Manhattan plots ───────────────────────────────────────────
echo ""
echo "[3] Generating Manhattan plots..."

Rscript $PIPELINE/_eQTL_manhattan.R \
    "$ANNOTATED_FILE" \
    "$LEAD_FILE" \
    "$OUTPUT_PREFIX" \
    "$FDR_THRESH"

echo "  Manhattan by SNP  : ${OUTPUT_PREFIX}.manhattan_by_SNP.pdf"
echo "  Manhattan by gene : ${OUTPUT_PREFIX}.manhattan_by_gene.pdf"

# ── Step 4: Boxplots ──────────────────────────────────────────────────
echo ""
echo "[4] Generating boxplots for top $TOP_N gene-SNP pairs..."

BOXPLOT_DIR=$OUTDIR/boxplots
mkdir -p $BOXPLOT_DIR

Rscript $PIPELINE/_eQTL_boxplot.R \
    "$TOP_PAIRS" \
    "$SNP_FILE" \
    "$EXPR_FILE" \
    "$BIM" \
    "$SNP_LOC" \
    "$BOXPLOT_DIR"

N_PLOTS=$(ls $BOXPLOT_DIR/*.pdf 2>/dev/null | wc -l)
echo "  Boxplots written : $N_PLOTS PDFs in $BOXPLOT_DIR"

# ── Final summary ─────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Matrix eQTL complete for : $TISSUE"
echo "  Results in               : $OUTDIR"
echo ""
echo "  Key output files:"
echo "  ${OUTPUT_PREFIX}.cis.txt                  (raw Matrix eQTL output)"
echo "  ${OUTPUT_PREFIX}.full_annotated.txt        (with chr/pos, flags)"
echo "  ${OUTPUT_PREFIX}.FDR${FDR_THRESH}.txt            (FDR-filtered)"
echo "  ${OUTPUT_PREFIX}.lead_snps.txt             (one lead SNP per gene)"
echo "  ${OUTPUT_PREFIX}.manhattan_by_SNP.pdf"
echo "  ${OUTPUT_PREFIX}.manhattan_by_gene.pdf"
echo "  $BOXPLOT_DIR/                              (top $TOP_N boxplots)"
echo ""
echo "  $(date)"
echo "============================================"
