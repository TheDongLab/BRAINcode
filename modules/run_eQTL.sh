#!/bin/bash
#SBATCH --job-name=run_eQTL_Combined
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/run_eQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/run_eQTL_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail
module load R

# ── Arguments ─────────────────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "ERROR: Missing tissue argument."
    echo "Usage: sbatch run_eQTL.sh \"Cerebellum\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  Matrix eQTL run for tissue : $TISSUE"
echo "  $(date)"
echo "============================================"

# ── Paths ─────────────────────────────────────────────────────────────
BASE=/home/zw529/donglab/data/target_ALS
PIPELINE=/home/zw529/donglab/pipelines/scripts/QTL
PLINK=$BASE/QTL/plink

# Working directory is the tissue's eQTL folder
INDIR=$BASE/$TISSUE_DIR/eQTL
OUTDIR=$INDIR/results
mkdir -p $OUTDIR

# Match the filenames from your 'ls' output
SNP_FILE=$INDIR/snp_${TISSUE_DIR}.txt
EXPR_FILE=$INDIR/expression_${TISSUE_DIR}.txt
COV_FILE=$INDIR/covariates_${TISSUE_DIR}_encoded.txt
GENE_LOC=$INDIR/gene_location.txt
SNP_LOC=$INDIR/snp_location.txt
BIM=$PLINK/joint_autosomes_filtered_bed.bim

# Output naming
OUTPUT_PREFIX=$OUTDIR/${TISSUE_DIR}_Combined_eQTL
CIS_FILE="${OUTPUT_PREFIX}.cis.txt"
FDR_THRESH=0.05
TOP_N=100

# ── Step 0: Pre-flight ────────────────────────────────────────────────
echo "[0] Verifying file existence and alignment..."
for f in "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" "$GENE_LOC" "$SNP_LOC"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File missing: $f"
        exit 1
    fi
done

# NEW: The Alignment Guard
# Counts columns in the three main matrices to ensure 100% agreement
S_N=$(head -n 1 "$SNP_FILE" | awk -F'\t' '{print NF-1}')
E_N=$(head -n 1 "$EXPR_FILE" | awk -F'\t' '{print NF-1}')
C_N=$(head -n 1 "$COV_FILE" | awk -F'\t' '{print NF-1}')

echo "Verification: SNPs($S_N), Expr($E_N), Covs($C_N)"

if [[ "$S_N" -ne "$E_N" ]] || [[ "$S_N" -ne "$C_N" ]]; then
    echo "FATAL ERROR: Matrices are not aligned! SNP:$S_N, EXPR:$E_N, COV:$C_N"
    exit 1
fi

# ── Step 1: Run Matrix eQTL ───────────────────────────────────────────
echo "[1] Running Matrix eQTL..."
Rscript $PIPELINE/_eQTL.R \
    "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" "$OUTPUT_PREFIX" "$GENE_LOC" "$SNP_LOC"

# ── Step 2: Post-processing ──────────────────────────────────────────
echo "[2] Post-processing..."
Rscript $PIPELINE/_eQTL_postprocess.R \
    "$CIS_FILE" "$SNP_LOC" "$GENE_LOC" "$OUTPUT_PREFIX" "$FDR_THRESH" "$TOP_N"

ANNOTATED_FILE="${OUTPUT_PREFIX}.full_annotated.txt"
LEAD_FILE="${OUTPUT_PREFIX}.lead_snps.txt"
TOP_PAIRS="${OUTPUT_PREFIX}.top_for_boxplot.txt"

# ── Step 3: Manhattan Plot ────────────────────────────────────────────
echo "[3] Generating Manhattan plot..."
Rscript $PIPELINE/_eQTL_manhattan.R \
    "$ANNOTATED_FILE" "$LEAD_FILE" "$OUTPUT_PREFIX" "$FDR_THRESH"

# ── Step 4: Boxplots ──────────────────────────────────────────────────
echo "[4] Generating boxplots for top $TOP_N pairs..."
# Passing OUTDIR directly so they print to results/ instead of a subfolder
Rscript $PIPELINE/_eQTL_boxplot.R \
    "$TOP_PAIRS" "$SNP_FILE" "$EXPR_FILE" "$BIM" "$SNP_LOC" "$OUTDIR"

# ── Final summary ─────────────────────────────────────────────────────
echo "============================================"
echo "  Run complete for : $TISSUE"
echo "  All outputs saved to: $OUTDIR"
echo "  $(date)"
echo "============================================"
