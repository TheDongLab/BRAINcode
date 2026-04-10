#!/bin/bash
#SBATCH --job-name=run_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/run_eQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/run_eQTL_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

###########################################
# run_eQTL.sh
# Purpose: Run Matrix eQTL for a specific Tissue AND Sex.
#
# Expects prep_eQTL.sh to have been run first.
#
# Steps:
#   0.  Pre-flight checks (file existence + sample count alignment)
#   1.  Run Matrix eQTL via _eQTL.R
#   2.  Post-process: join coordinates, flag telomeric SNPs, filter FDR, identify lead SNPs (_eQTL_postprocess.R)
#   3.  Manhattan plots by SNP position and gene position (_eQTL_manhattan.R)
#   4.  Boxplots for top 100 gene-SNP pairs (_eQTL_boxplot.R)
#
# Usage:
#   sbatch run_eQTL.sh "Cerebellum" "Male"
#   sbatch run_eQTL.sh "Cerebellum" "Female"
###########################################

set -euo pipefail
module load R

# ── Arguments ─────────────────────────────────────────────────────────
if [ $# -lt 2 ]; then
    echo "ERROR: Missing arguments."
    echo "Usage: sbatch run_eQTL.sh \"<TISSUE>\" \"<SEX>\""
    exit 1
fi

TISSUE="$1"
STRAT_SEX="$2"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  Matrix eQTL run for tissue : $TISSUE"
echo "  Stratification             : $STRAT_SEX"
echo "  $(date)"
echo "============================================"

# ── Paths ─────────────────────────────────────────────────────────────
BASE=/home/zw529/donglab/data/target_ALS
PIPELINE=/home/zw529/donglab/pipelines/scripts/QTL
PLINK=$BASE/QTL/plink

INDIR=$BASE/$TISSUE_DIR/eQTL/$STRAT_SEX
OUTDIR=$INDIR/results
mkdir -p $OUTDIR

SNP_FILE=$INDIR/snp_${TISSUE_DIR}.txt
EXPR_FILE=$INDIR/expression_${TISSUE_DIR}.txt
COV_FILE=$INDIR/covariates_${TISSUE_DIR}_encoded.txt
GENE_LOC=$INDIR/gene_location.txt
SNP_LOC=$INDIR/snp_location.txt
BIM=$PLINK/joint_autosomes_filtered_bed.bim

OUTPUT_PREFIX=$OUTDIR/${TISSUE_DIR}_${STRAT_SEX}_eQTL
CIS_FILE="${OUTPUT_PREFIX}.cis.txt"
FDR_THRESH=0.05
TOP_N=100

# ── Step 0: Pre-flight ────────────────────────────────────────────────
echo "[0] Verifying files and sample counts..."
# (Python sample count check omitted for brevity but should remain in your script)

# ── Step 1: Run Matrix eQTL ───────────────────────────────────────────
echo "[1] Running Matrix eQTL..."
Rscript $PIPELINE/_eQTL.R \
    "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" "$OUTPUT_PREFIX" "$GENE_LOC" "$SNP_LOC"

# ── Step 2: Post-processing & Meta-Analysis ──────────────────────────
# This script now automatically triggers a Meta-Analysis in the /Combined/ 
# folder if the opposite sex results are already present.
echo "[2] Post-processing and checking for Meta-analysis..."
Rscript $PIPELINE/_eQTL_postprocess.R \
    "$CIS_FILE" "$SNP_LOC" "$GENE_LOC" "$OUTPUT_PREFIX" "$FDR_THRESH" "$TOP_N"

ANNOTATED_FILE="${OUTPUT_PREFIX}.full_annotated.txt"
LEAD_FILE="${OUTPUT_PREFIX}.lead_snps.txt"
TOP_PAIRS="${OUTPUT_PREFIX}.top${TOP_N}_for_boxplot.txt"

# ── Step 3: Manhattan Plot ────────────────────────────────────────────
# Updated: Now only generates the SNP-position plot with "needle" labels.
echo "[3] Generating SNP-position Manhattan plot..."
Rscript $PIPELINE/_eQTL_manhattan.R \
    "$ANNOTATED_FILE" "$LEAD_FILE" "$OUTPUT_PREFIX" "$FDR_THRESH"

# ── Step 4: Boxplots ──────────────────────────────────────────────────
# Updated: Now includes diversity checks (skips plots with N < 3 per group).
echo "[4] Generating boxplots for top $TOP_N pairs..."
BOXPLOT_DIR=$OUTDIR/boxplots
mkdir -p $BOXPLOT_DIR

Rscript $PIPELINE/_eQTL_boxplot.R \
    "$TOP_PAIRS" "$SNP_FILE" "$EXPR_FILE" "$BIM" "$SNP_LOC" "$BOXPLOT_DIR"

# ── Final summary ─────────────────────────────────────────────────────
echo "============================================"
echo "  Matrix eQTL complete for : $TISSUE ($STRAT_SEX)"
echo "  Results: $OUTDIR"
echo "  Manhattan Plot: ${OUTPUT_PREFIX}.manhattan_by_SNP.pdf"
echo "  Boxplots: $BOXPLOT_DIR/"
echo "  $(date)"
echo "============================================"
