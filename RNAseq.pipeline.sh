#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --cpus-per-task=4
#SBATCH --mem=99G
#SBATCH --time=24:00:00
#SBATCH -p day
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# ── Meta RNAseq submission script ─────────────────────────────────────────────
# Discovers all RNAseq samples across target experiments and runs _RNAseq.sh
# on each sample as a SLURM array job.
#
# For most experiments: discovers samples from RNAseq/Raw/ subdirectories
#                       and writes output to RNAseq/Processed/.
# For target_ALS:       discovers samples from RNAseq/Processed/ subdirectories
#                       (BAM already present; steps 1-4 are skipped via status
#                       files, and the BAM is symlinked to the expected name).
#
# Lives in: $HOME/donglab/pipelines/scripts/rnaseq/
#
# Usage:
#   bash RNAseq.pipeline.sh            # count samples and self-submit
#   bash RNAseq.pipeline.sh --dry-run  # preview without submitting
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR="/home/zw529/donglab/data"
PIPELINE_DIR="$HOME/donglab/pipelines/scripts/rnaseq"
SCRIPT="$PIPELINE_DIR/_RNAseq.sh"

# ── Target experiments ────────────────────────────────────────────────────────
# Format: "EXPERIMENT_NAME|raw"       --> discover from RNAseq/Raw/, output to RNAseq/Processed/
#         "EXPERIMENT_NAME|processed" --> discover from RNAseq/Processed/ (BAM already present)
EXPERIMENTS=(
    "AMPALS_GSE124439|raw"
    "AMPALS_GSE219281|raw"
    "GEO_GSE153960|raw"
    "target_ALS|processed"
)

# ── Shared function: build ordered master sample list into MASTER array ────────
# Each entry: "EXP|tissue|sample|mode"
build_master() {
    MASTER=()
    for entry in "${EXPERIMENTS[@]}"; do
        EXP=$(echo "$entry"  | cut -d'|' -f1)
        MODE=$(echo "$entry" | cut -d'|' -f2)
        EXP_DIR="$DATA_DIR/$EXP"
        [[ -d "$EXP_DIR" ]] || continue

        if [[ "$MODE" == "raw" ]]; then
            # Discover sample subdirectories inside RNAseq/Raw/
            while IFS= read -r raw_dir; do
                rnaseq_dir=$(dirname "$raw_dir")
                tissue=$(basename "$(dirname "$rnaseq_dir")")
                for sample_dir in "$raw_dir"/*/; do
                    [[ -d "$sample_dir" ]] || continue
                    sample=$(basename "$sample_dir")
                    MASTER+=("${EXP}|${tissue}|${sample}|raw")
                done
            done < <(find "$EXP_DIR" -type d -path '*/RNAseq/Raw' 2>/dev/null | sort)

        elif [[ "$MODE" == "processed" ]]; then
            # Discover sample subdirectories inside RNAseq/Processed/
            while IFS= read -r proc_dir; do
                rnaseq_dir=$(dirname "$proc_dir")
                tissue=$(basename "$(dirname "$rnaseq_dir")")
                for sample_dir in "$proc_dir"/*/; do
                    [[ -d "$sample_dir" ]] || continue
                    sample=$(basename "$sample_dir")
                    MASTER+=("${EXP}|${tissue}|${sample}|processed")
                done
            done < <(find "$EXP_DIR" -type d -path '*/RNAseq/Processed' 2>/dev/null | sort)
        fi
    done
}

# ── Array task: process one sample ────────────────────────────────────────────
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then

    build_master

    ENTRY="${MASTER[$((SLURM_ARRAY_TASK_ID - 1))]}"
    EXP=$(echo "$ENTRY"    | cut -d'|' -f1)
    TISSUE=$(echo "$ENTRY" | cut -d'|' -f2)
    SAMPLE=$(echo "$ENTRY" | cut -d'|' -f3)
    MODE=$(echo "$ENTRY"   | cut -d'|' -f4)

    echo "=============================="
    echo "Task   : $SLURM_ARRAY_TASK_ID"
    echo "Exp    : $EXP"
    echo "Tissue : $TISSUE"
    echo "Sample : $SAMPLE"
    echo "Mode   : $MODE"

    if [[ "$MODE" == "raw" ]]; then

        RAW_DIR="$DATA_DIR/$EXP/$TISSUE/RNAseq/Raw/$SAMPLE"
        OUTDIR="$DATA_DIR/$EXP/$TISSUE/RNAseq/Processed/$SAMPLE"

        FASTQ1=$(ls "$RAW_DIR"/*_1.fastq.gz 2>/dev/null | head -1 || true)
        FASTQ2=$(ls "$RAW_DIR"/*_2.fastq.gz 2>/dev/null | head -1 || true)

        if [[ -z "$FASTQ1" || -z "$FASTQ2" ]]; then
            echo "ERROR: No paired FASTQs found in $RAW_DIR — skipping."
            exit 1
        fi

        echo "FASTQ1 : $FASTQ1"
        echo "FASTQ2 : $FASTQ2"
        echo "Output : $OUTDIR"
        echo "=============================="

        mkdir -p "$OUTDIR"
        cd "$PIPELINE_DIR"
        bash "$SCRIPT" "$FASTQ1" "$FASTQ2" "$OUTDIR" \
            > "$OUTDIR/output.out" \
            2> "$OUTDIR/output.err"

    elif [[ "$MODE" == "processed" ]]; then

        OUTDIR="$DATA_DIR/$EXP/$TISSUE/RNAseq/Processed/$SAMPLE"
        EXPECTED_BAM="$OUTDIR/STAR.Aligned.sortedByCoord.out.bam"

        echo "Output : $OUTDIR"
        echo "=============================="

        # Find the existing BAM (*.sorted.bam or *.final.bam)
        EXISTING_BAM=$(ls "$OUTDIR"/*.sorted.bam "$OUTDIR"/*.final.bam 2>/dev/null | head -1 || true)
        if [[ -z "$EXISTING_BAM" ]]; then
            echo "ERROR: No .sorted.bam or .final.bam found in $OUTDIR — skipping."
            exit 1
        fi

        # Symlink BAM and BAI to the names _RNAseq.sh expects
        if [[ ! -f "$EXPECTED_BAM" ]]; then
            ln -s "$EXISTING_BAM" "$EXPECTED_BAM"
            echo "Symlinked: $(basename "$EXISTING_BAM") -> STAR.Aligned.sortedByCoord.out.bam"
        fi
        if [[ ! -f "${EXPECTED_BAM}.bai" ]]; then
            EXISTING_BAI="${EXISTING_BAM}.bai"
            if [[ -f "$EXISTING_BAI" ]]; then
                ln -s "$EXISTING_BAI" "${EXPECTED_BAM}.bai"
                echo "Symlinked: $(basename "$EXISTING_BAI") -> STAR.Aligned.sortedByCoord.out.bam.bai"
            else
                echo "WARNING: No .bai found — samtools index will run in step 5."
            fi
        fi

        # Pre-touch status files for steps 1-4 so _RNAseq.sh skips them
        for step in fastqc trim kpal mapping; do
            touch "$OUTDIR/.status.RNAseq.${step}"
        done
        echo "Pre-touched status files for steps 1-4 (fastqc, trim, kpal, mapping)."

        # Call _RNAseq.sh with dummy FASTQ args (steps 1-4 will be skipped)
        cd "$PIPELINE_DIR"
        bash "$SCRIPT" "SKIP" "SKIP" "$OUTDIR" \
            > "$OUTDIR/output.out" \
            2> "$OUTDIR/output.err"
    fi

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
