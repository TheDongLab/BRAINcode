#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH -p day

# ── Meta RNAseq submission script ─────────────────────────────────────────────
# Discovers all sample_list.txt files across target experiments and runs
# _RNAseq.sh on each sample as a SLURM array job.
#
# Lives in: $HOME/donglab/pipelines/scripts/rnaseq/
#
# Usage:
#   bash run_RNAseq_all.sh            # count samples and self-submit
#
# Dry-run (preview samples without submitting):
#   bash run_RNAseq_all.sh --dry-run
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR="/home/zw529/donglab/data"
PIPELINE_DIR="$HOME/donglab/pipelines/scripts/rnaseq"
SCRIPT="$PIPELINE_DIR/_RNAseq.sh"

# ── Target experiments ────────────────────────────────────────────────────────
EXPERIMENTS=(
    "AMPALS_GSE124439"
    "AMPALS_GSE219281"
    "GEO_GSE153960"
    "target_ALS"
)

# ── Shared function: build ordered master sample list into MASTER array ────────
build_master() {
    MASTER=()
    for EXP in "${EXPERIMENTS[@]}"; do
        EXP_DIR="$DATA_DIR/$EXP"
        [[ -d "$EXP_DIR" ]] || continue
        while IFS= read -r sample_list; do
            rnaseq_dir=$(dirname "$sample_list")
            tissue=$(basename "$(dirname "$rnaseq_dir")")
            while IFS= read -r sample; do
                [[ -z "$sample" ]] && continue
                MASTER+=("${EXP}|${tissue}|${sample}")
            done < "$sample_list"
        done < <(find "$EXP_DIR" -type f -path '*/RNAseq/sample_list.txt' 2>/dev/null | sort)
    done
}

# ── Array task: process one sample ────────────────────────────────────────────
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then

    build_master

    ENTRY="${MASTER[$((SLURM_ARRAY_TASK_ID - 1))]}"
    EXP=$(echo "$ENTRY"    | cut -d'|' -f1)
    TISSUE=$(echo "$ENTRY" | cut -d'|' -f2)
    SAMPLE=$(echo "$ENTRY" | cut -d'|' -f3)

    RAW_DIR="$DATA_DIR/$EXP/$TISSUE/RNAseq/Raw/$SAMPLE"
    OUTDIR="$DATA_DIR/$EXP/$TISSUE/RNAseq/processed/$SAMPLE"

    FASTQ1=$(ls "$RAW_DIR"/*_1.fastq.gz 2>/dev/null | head -1 || true)
    FASTQ2=$(ls "$RAW_DIR"/*_2.fastq.gz 2>/dev/null | head -1 || true)

    if [[ -z "$FASTQ1" || -z "$FASTQ2" ]]; then
        echo "ERROR: No paired FASTQs found in $RAW_DIR — skipping."
        exit 1
    fi

    mkdir -p "$OUTDIR"

    echo "=============================="
    echo "Task   : $SLURM_ARRAY_TASK_ID"
    echo "Exp    : $EXP"
    echo "Tissue : $TISSUE"
    echo "Sample : $SAMPLE"
    echo "FASTQ1 : $FASTQ1"
    echo "FASTQ2 : $FASTQ2"
    echo "Output : $OUTDIR"
    echo "=============================="

    cd "$PIPELINE_DIR"
    bash "$SCRIPT" "$FASTQ1" "$FASTQ2" "$OUTDIR" \
        > "$OUTDIR/output.out" \
        2> "$OUTDIR/output.err"

    exit 0
fi

# ── Initial run: count samples and self-submit as array ───────────────────────
echo "Scanning for samples..."
build_master
COUNT=${#MASTER[@]}
echo "Found $COUNT samples."

if [[ "$COUNT" -eq 0 ]]; then
    echo "No samples found. Exiting."
    exit 1
fi

if $DRY_RUN; then
    for i in "${!MASTER[@]}"; do
        echo "  [$((i+1))] ${MASTER[$i]}"
    done
    echo ""
    echo "Would submit: --array=1-${COUNT}%50"
    exit 0
fi

echo "Submitting array job: 1-${COUNT} (max 50 concurrent)..."
sbatch --array="1-${COUNT}%50" "$0"
echo "Done. Monitor with: squeue -u $USER"
