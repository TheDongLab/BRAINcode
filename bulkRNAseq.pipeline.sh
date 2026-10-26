#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --cpus-per-task=4
#SBATCH --mem=99G
#SBATCH --time=24:00:00
#SBATCH -p day
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Submit one RNA-seq worker per sample.

set -euo pipefail

###########################################
# Arguments
###########################################
DRY_RUN=false
TARGET_USER=""

while (($#)); do
    case "$1" in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        *)
            TARGET_USER="$1"
            shift
            ;;
    esac
done

TARGET_USER="${TARGET_USER:-$USER}"

###########################################
# Paths and settings
###########################################
DATA_DIR="/home/${TARGET_USER}/donglab/data"
PIPELINE_DIR="/home/${TARGET_USER}/donglab/pipelines/scripts/rnaseq"
SCRIPT="$PIPELINE_DIR/_RNAseq.sh"

ANSWERALS_TISSUE="iPSC_Motor_Neuron"
ANSWERALS_RNASEQ_DIR="$DATA_DIR/answer_ALS/${ANSWERALS_TISSUE}/RNAseq"
ANSWERALS_PROCESSED_DIR="$ANSWERALS_RNASEQ_DIR/Processed"

###########################################
# Experiments to include
#
# Keep TargetALS as the current default.
# Uncomment AnswerALS when ready to include that cohort.
###########################################
EXPERIMENTS=(
    "target_ALS|processed"
    # "answer_ALS|answerals"
)

###########################################
# Build deterministic ordered sample list
# Entries are:
#   EXP|TISSUE|SAMPLE|MODE|BAM
###########################################
build_master() {
    MASTER=()

    for entry in "${EXPERIMENTS[@]}"; do
        EXP="${entry%%|*}"
        REST="${entry#*|}"
        MODE="$REST"

        if [[ "$EXP" == "answer_ALS" && "$MODE" == "answerals" ]]; then
            [[ -d "$ANSWERALS_PROCESSED_DIR" ]] || {
                echo "[ERROR] AnswerALS RNA-seq directory not found: $ANSWERALS_PROCESSED_DIR" >&2
                return 1
            }

            while IFS= read -r donor_dir; do
                [[ -d "$donor_dir" ]] || continue
                donor="$(basename "$donor_dir")"

                # This control was excluded from the agreed 1,101-BAM cohort.
                [[ "$donor" == "CTRL-NEUEU392AE8" ]] && continue

                # Select exactly one original BAM. Exclude release5, remap,
                # and STAR-generated BAMs; fail on ambiguous donors.
                matches=()
                while IFS= read -r bam; do
                    [[ -f "$bam" ]] || continue
                    case "$(basename "$bam")" in
                        *release5*|*remap*|STAR.*) continue ;;
                    esac
                    matches+=("$bam")
                done < <(find "$donor_dir" -maxdepth 1 -type f -name "*.sorted.bam" | sort)

                if ((${#matches[@]} != 1)); then
                    echo "[ERROR] Expected exactly one selected original BAM for $donor; found ${#matches[@]}." >&2
                    printf '  %s\n' "${matches[@]:-}" >&2
                    return 1
                fi

                MASTER+=("answer_ALS|${ANSWERALS_TISSUE}|${donor}|answerals|${matches[0]}")
            done < <(find "$ANSWERALS_PROCESSED_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

        elif [[ "$MODE" == "raw" ]]; then
            EXP_DIR="$DATA_DIR/$EXP"
            [[ -d "$EXP_DIR" ]] || continue

            while IFS= read -r raw_dir; do
                rnaseq_dir="$(dirname "$raw_dir")"
                tissue="$(basename "$(dirname "$rnaseq_dir")")"

                for sample_dir in "$raw_dir"/*/; do
                    [[ -d "$sample_dir" ]] || continue
                    sample="$(basename "$sample_dir")"
                    MASTER+=("${EXP}|${tissue}|${sample}|raw|")
                done
            done < <(find "$EXP_DIR" -type d -path '*/RNAseq/Raw' 2>/dev/null | sort)

        elif [[ "$MODE" == "processed" ]]; then
            EXP_DIR="$DATA_DIR/$EXP"
            [[ -d "$EXP_DIR" ]] || continue

            while IFS= read -r proc_dir; do
                rnaseq_dir="$(dirname "$proc_dir")"
                tissue="$(basename "$(dirname "$rnaseq_dir")")"

                for sample_dir in "$proc_dir"/*/; do
                    [[ -d "$sample_dir" ]] || continue
                    sample="$(basename "$sample_dir")"
                    MASTER+=("${EXP}|${tissue}|${sample}|processed|")
                done
            done < <(find "$EXP_DIR" -type d -path '*/RNAseq/Processed' 2>/dev/null | sort)
        fi
    done
}

###########################################
# Array task: process one sample
###########################################
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    build_master

    TASK_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
    if ((TASK_INDEX < 0 || TASK_INDEX >= ${#MASTER[@]})); then
        echo "[ERROR] Array task ID $SLURM_ARRAY_TASK_ID is outside the sample list." >&2
        exit 1
    fi

    ENTRY="${MASTER[$TASK_INDEX]}"
    IFS='|' read -r EXP TISSUE SAMPLE MODE ANSWERALS_BAM <<< "$ENTRY"

    echo "=============================="
    echo "Task   : $SLURM_ARRAY_TASK_ID"
    echo "User   : $TARGET_USER"
    echo "Exp    : $EXP"
    echo "Tissue : $TISSUE"
    echo "Sample : $SAMPLE"
    echo "Mode   : $MODE"

    if [[ "$MODE" == "answerals" ]]; then
        OUTDIR="$ANSWERALS_PROCESSED_DIR/$SAMPLE"

        echo "BAM    : $ANSWERALS_BAM"
        echo "Output : $OUTDIR"
        echo "=============================="

        mkdir -p "$OUTDIR"
        cd "$PIPELINE_DIR"

        bash "$SCRIPT" "SKIP" "SKIP" "$OUTDIR" "$TARGET_USER" "$ANSWERALS_BAM" \
            > "$OUTDIR/output.out" \
            2> "$OUTDIR/output.err"

    elif [[ "$MODE" == "raw" ]]; then
        RAW_DIR="$DATA_DIR/$EXP/$TISSUE/RNAseq/Raw/$SAMPLE"
        OUTDIR="$DATA_DIR/$EXP/$TISSUE/RNAseq/Processed/$SAMPLE"

        FASTQ1=$(ls "$RAW_DIR"/*_1.fastq.gz 2>/dev/null | head -1 || true)
        FASTQ2=$(ls "$RAW_DIR"/*_2.fastq.gz 2>/dev/null | head -1 || true)

        if [[ -z "$FASTQ1" || -z "$FASTQ2" ]]; then
            echo "[ERROR] No paired FASTQs found in $RAW_DIR." >&2
            exit 1
        fi

        echo "FASTQ1 : $FASTQ1"
        echo "FASTQ2 : $FASTQ2"
        echo "Output : $OUTDIR"
        echo "=============================="

        mkdir -p "$OUTDIR"
        cd "$PIPELINE_DIR"

        bash "$SCRIPT" "$FASTQ1" "$FASTQ2" "$OUTDIR" "$TARGET_USER" \
            > "$OUTDIR/output.out" \
            2> "$OUTDIR/output.err"

    elif [[ "$MODE" == "processed" ]]; then
        OUTDIR="$DATA_DIR/$EXP/$TISSUE/RNAseq/Processed/$SAMPLE"
        EXPECTED_BAM="$OUTDIR/STAR.Aligned.sortedByCoord.out.bam"

        echo "Output : $OUTDIR"
        echo "=============================="

        # Preserve the existing TargetALS BAM discovery and call behavior.
        EXISTING_BAM=""
        while IFS= read -r candidate; do
            if [[ -f "$candidate" && ! "$candidate" =~ remap ]]; then
                EXISTING_BAM="$candidate"
                break
            fi
        done < <(find "$OUTDIR" -maxdepth 1 -type f -name "*.bam" 2>/dev/null)

        if [[ -z "$EXISTING_BAM" ]]; then
            echo "[ERROR] No valid, non-remapped BAM found in $OUTDIR." >&2
            exit 1
        fi

        ln -sf "$EXISTING_BAM" "$EXPECTED_BAM"
        echo "Symlinked: $(basename "$EXISTING_BAM") -> STAR.Aligned.sortedByCoord.out.bam"

        EXISTING_BAI=""
        if [[ -f "${EXISTING_BAM}.bai" ]]; then
            EXISTING_BAI="${EXISTING_BAM}.bai"
        elif [[ -f "${EXISTING_BAM%.bam}.bai" ]]; then
            EXISTING_BAI="${EXISTING_BAM%.bam}.bai"
        fi

        if [[ -n "$EXISTING_BAI" ]]; then
            ln -sf "$EXISTING_BAI" "${EXPECTED_BAM}.bai"
            echo "Symlinked: $(basename "$EXISTING_BAI") -> STAR.Aligned.sortedByCoord.out.bam.bai"
        else
            echo "WARNING: No .bai found — samtools index will run when needed."
        fi

        for step in fastqc trim kpal mapping; do
            touch "$OUTDIR/.status.RNAseq.${step}"
        done

        cd "$PIPELINE_DIR"
        bash "$SCRIPT" "SKIP" "SKIP" "$OUTDIR" "$TARGET_USER" \
            > "$OUTDIR/output.out" \
            2> "$OUTDIR/output.err"
    fi

    exit 0
fi

###########################################
# Initial run: build list, preview or submit
###########################################
echo "Scanning data directory for user: $TARGET_USER"
build_master

COUNT=${#MASTER[@]}
echo "Found $COUNT samples."

if [[ "$COUNT" -eq 0 ]]; then
    echo "[ERROR] No samples found. Check EXPERIMENTS and the data paths." >&2
    exit 1
fi

if "$DRY_RUN"; then
    for i in "${!MASTER[@]}"; do
        echo "  [$((i + 1))] ${MASTER[$i]}"
    done
    echo
    echo "Would submit array: --array=1-${COUNT}%99"
    exit 0
fi

echo "Submitting array job: 1-${COUNT} (max 99 concurrent)..."
ARRAY_SUBMISSION=$(sbatch --parsable --array="1-${COUNT}%99" "$0" "$TARGET_USER")
ARRAY_JOB_ID="${ARRAY_SUBMISSION%%;*}"
echo "Submitted array job $ARRAY_JOB_ID."

echo "Monitor with: squeue -u $TARGET_USER"
