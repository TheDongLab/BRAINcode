#!/bin/bash
#SBATCH --job-name=run_smr_heidi_plot
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/run_smr_heidi_plot.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/run_smr_heidi_plot.err
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

TYPE="${1:-}"; ROOT="/home/zw529/donglab/data/target_ALS"; MR_DIR="${ROOT}/MR"
METADATA="${ROOT}/targetALS_rnaseq_metadata.csv"
GTF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf"
RAW_FILE="${ROOT}/QTL/plink/joint_all_chrs_matrixEQTL.raw"
SMR_FDR_THRESHOLD=0.05; FLANK_FRAC=.20; FLANK_MIN=2000; FLANK_MAX=20000; BIGWIG_BIN_BP=25; MAX_TRACK_POINTS=2500
CIRC_COORD_TOL=3; MIN_MEAN_JUNCTION_RPM=.01; MAX_JUNCTIONS=30; MAX_HITS=0; TISSUE_FILTER=""; OVERWRITE=0; DPI=180

set -euo pipefail
case "$TYPE" in eQTL|sQTL|cQTL) ;; *) echo "Usage: sbatch $0 {eQTL|sQTL|cQTL}"; exit 1 ;; esac

OUTDIR="${MR_DIR}/SMR_HEIDI_raw_plots/${TYPE}"; mkdir -p "$OUTDIR"
for f in "$METADATA" "$GTF" "$RAW_FILE"; do [[ -f "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done

module load deepTools 2>/dev/null || { echo "ERROR: could not load deepTools"; exit 1; }
module load poppler/25.07.0-GCC-13.3.0 2>/dev/null || { echo "ERROR: could not load Poppler"; exit 1; }

PYTHON="$(command -v python)"
for x in pdfseparate pdftotext samtools; do command -v "$x" >/dev/null || { echo "ERROR: $x unavailable"; exit 1; }; done
"$PYTHON" -c 'import numpy,matplotlib,pyBigWig' || { echo "ERROR: Python plotting stack unavailable"; exit 1; }

export TYPE ROOT MR_DIR METADATA OUTDIR GTF RAW_FILE SMR_FDR_THRESHOLD FLANK_FRAC FLANK_MIN FLANK_MAX BIGWIG_BIN_BP MAX_TRACK_POINTS
export CIRC_COORD_TOL MIN_MEAN_JUNCTION_RPM MAX_JUNCTIONS MAX_HITS TISSUE_FILTER OVERWRITE DPI

SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
"$PYTHON" "$SCRIPT_DIR/smr_heidi_plot.py"
