#!/bin/bash
#SBATCH --job-name=run_sQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/run_sQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/run_sQTL_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail
module load R

# ── Arguments ─────────────────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "ERROR: Missing tissue argument."
    echo "Usage: sbatch run_sQTL.sh \"Cerebellum\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  Matrix sQTL run for tissue : $TISSUE"
echo "  $(date)"
echo "============================================"

# ── Paths ─────────────────────────────────────────────────────────────
BASE=/home/zw529/donglab/data/target_ALS
PIPELINE=/home/zw529/donglab/pipelines/scripts/QTL
PLINK=$BASE/QTL/plink

# Working directory is the tissue's sQTL folder
INDIR=$BASE/$TISSUE_DIR/sQTL
OUTDIR=$INDIR/results
mkdir -p $OUTDIR

# Match the file names from your successful prep_sQTL outputs
SNP_FILE=$INDIR/snp_${TISSUE_DIR}.txt
SPLICING_FILE=$INDIR/splicing_${TISSUE_DIR}.txt
COV_FILE=$INDIR/covariates_${TISSUE_DIR}_encoded.txt
SPLICING_LOC=$INDIR/splicing_location.txt
SNP_LOC=$INDIR/snp_location.txt
BIM=$PLINK/joint_all_chrs_filtered_bed.bim

# Output naming (Keeping files flat in $OUTDIR)
OUTPUT_PREFIX=$OUTDIR/${TISSUE_DIR}_sQTL
CIS_FILE="${OUTPUT_PREFIX}.cis.txt"
FDR_THRESH=0.05
TOP_N=1000000

# ── Step 0: Pre-flight & Location Harmonization ──────────────────────
echo "[0] Verifying file existence and alignment..."
for f in "$SNP_FILE" "$SPLICING_FILE" "$COV_FILE" "$SPLICING_LOC" "$SNP_LOC"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File missing: $f"
        exit 1
    fi
done

# Harmonization Guard: Clean chromosome suffixes from splicing locations
# This strips trailing strandedness signs (e.g., 'chr10:+' or 'chr10:-' -> 'chr10')
# so the chromosomes perfectly align with the structural map in snp_location.txt
echo "[0.1] Standardizing chromosome labels in splicing location file..."
SPLICING_LOC_BAK="${SPLICING_LOC}.bak"

# Only make a fresh backup if it doesn't already exist to preserve the pristine original raw files
if [ ! -f "$SPLICING_LOC_BAK" ]; then
    cp "$SPLICING_LOC" "$SPLICING_LOC_BAK"
fi

# Clean column 2 while preserving tab structure safely
sed -E 's/(chr[0-9XY]+):[\+\-]/\1/g' "$SPLICING_LOC_BAK" > "$SPLICING_LOC"
echo "      -> Chromosome strings cleaned in $SPLICING_LOC"

# Alignment Guard
S_N=$(head -n 1 "$SNP_FILE" | awk -F'\t' '{print NF-1}')
SPL_N=$(head -n 1 "$SPLICING_FILE" | awk -F'\t' '{print NF-1}')
C_N=$(head -n 1 "$COV_FILE" | awk -F'\t' '{print NF-1}')

echo "Verification: SNPs($S_N), Splicing($SPL_N), Covs($C_N)"

if [[ "$S_N" -ne "$SPL_N" ]] || [[ "$S_N" -ne "$C_N" ]]; then
    echo "FATAL ERROR: Matrices are not aligned! SNP:$S_N, SPLICING:$SPL_N, COV:$C_N"
    exit 1
fi

# ── Step 1: Run Matrix sQTL ───────────────────────────────────────────
echo "[1] Running Matrix sQTL..."
Rscript $PIPELINE/_sQTL.R \
    "$SNP_FILE" "$SPLICING_FILE" "$COV_FILE" "$OUTPUT_PREFIX" "$SPLICING_LOC" "$SNP_LOC"

# ── Step 2: Post-processing ──────────────────────────────────────────
echo "[2] Post-processing..."
Rscript $PIPELINE/_sQTL_postprocess.R \
    "$CIS_FILE" "$SNP_LOC" "$SPLICING_LOC" "$OUTPUT_PREFIX" "$FDR_THRESH" "$TOP_N"

ANNOTATED_FILE="${OUTPUT_PREFIX}.full_annotated.txt"
LEAD_FILE="${OUTPUT_PREFIX}.lead_snps.txt"
TOP_PAIRS="${OUTPUT_PREFIX}.top_for_boxplot.txt"

# ── Step 3: Manhattan Plot ────────────────────────────────────────────
echo "[3] Generating Manhattan plot..."
Rscript $PIPELINE/_sQTL_manhattan.R \
    "$ANNOTATED_FILE" "$LEAD_FILE" "$OUTPUT_PREFIX" "$FDR_THRESH"

# ── Step 4: Boxplots ──────────────────────────────────────────────────
echo "[4] Generating boxplots for all sig. SNPs..."
Rscript $PIPELINE/_sQTL_boxplot.R \
    "$TOP_PAIRS" "$SNP_FILE" "$SPLICING_FILE" "$COV_FILE" "$SNP_LOC" "$OUTDIR" "$TISSUE_DIR"

# ── Step 5: Cleanup Directory Sprawl ──────────────────────────────────
# If a folder was created with the prefix name, move contents up and delete it
if [ -d "${OUTPUT_PREFIX}" ]; then
    echo "[5] Cleaning up redundant subdirectories..."
    mv "${OUTPUT_PREFIX}"/* "$OUTDIR/" 2>/dev/null || true
    rmdir "${OUTPUT_PREFIX}" 2>/dev/null || true
fi

# ── Final summary ─────────────────────────────────────────────────────
echo "============================================"
echo "  Run complete for : $TISSUE"
echo "  All outputs saved to: $OUTDIR"
echo "  $(date)"
echo "============================================"
