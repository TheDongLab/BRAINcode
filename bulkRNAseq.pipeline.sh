#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --cpus-per-task=4
#SBATCH --mem=99G
#SBATCH --time=24:00:00
#SBATCH -p day
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Submit one RNA-seq worker per selected sample.

set -euo pipefail

###########################################
# Arguments
###########################################
DRY_RUN=false
TARGET_USER=""
SELECTION_MANIFEST=""

while (($#)); do
    case "$1" in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --manifest)
            [[ $# -ge 2 ]] || {
                echo "[ERROR] --manifest requires a file path." >&2
                exit 1
            }
            SELECTION_MANIFEST="$2"
            shift 2
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

GENOME_BASE="/home/${TARGET_USER}/donglab/references/genome/Homo_sapiens/UCSC/hg38"
ANSWERALS_CHROMINFO="$GENOME_BASE/Sequence/STAR/genome.sizes.nochr"

QC_BED=""

cleanup() {
    [[ -z "$QC_BED" ]] || rm -f "$QC_BED"
}
trap cleanup EXIT

###########################################
# Experiments to include
#
# Keep AnswerALS as the current selection.
# Uncomment TargetALS when you want to include that cohort.
###########################################
EXPERIMENTS=(
    # "target_ALS|processed"
    "answer_ALS|answerals"
)

###########################################
# QC helpers used only if multiple candidates
# have matching chimeric remaps.
###########################################
ensure_qc_tools() {
    if ! command -v samtools >/dev/null 2>&1; then
        set +u
        source "/home/${TARGET_USER}/donglab/pipelines/modules/miniconda3/etc/profile.d/conda.sh"
        conda activate RNAseq
        set -u
    fi

    command -v samtools >/dev/null 2>&1 || {
        echo "[ERROR] samtools is unavailable; cannot apply the BAM QC tiebreaker." >&2
        return 1
    }

    [[ -s "$ANSWERALS_CHROMINFO" ]] || {
        echo "[ERROR] AnswerALS canonical chromosome sizes are missing or empty: $ANSWERALS_CHROMINFO" >&2
        return 1
    }

    if [[ -z "$QC_BED" ]]; then
        QC_BED="$(mktemp "${TMPDIR:-/tmp}/answerALS_selection_XXXXXX.bed")"
        awk 'BEGIN { FS="[ \t]+"; OFS="\t" } NF >= 2 { print $1, 0, $2 }' \
            "$ANSWERALS_CHROMINFO" > "$QC_BED"
    fi
}

###########################################
# Build or load deterministic ordered sample list
#
# AnswerALS entries are:
#   EXP|TISSUE|DONOR|MODE|BAM|CHIMERIC_REMAP|QC_COUNT|SELECTION_REASON
###########################################
build_master() {
    MASTER=()

    # The submitted array reads the exact sample choices made before submission.
    if [[ -n "$SELECTION_MANIFEST" ]]; then
        [[ -r "$SELECTION_MANIFEST" ]] || {
            echo "[ERROR] Selection manifest is not readable: $SELECTION_MANIFEST" >&2
            return 1
        }
        mapfile -t MASTER < "$SELECTION_MANIFEST"
        ((${#MASTER[@]} > 0)) || {
            echo "[ERROR] Selection manifest is empty: $SELECTION_MANIFEST" >&2
            return 1
        }
        return 0
    fi

    for entry in "${EXPERIMENTS[@]}"; do
        EXP="${entry%%|*}"
        REST="${entry#*|}"
        MODE="$REST"

        if [[ "$EXP" == "answer_ALS" && "$MODE" == "answerals" ]]; then
            [[ -d "$ANSWERALS_PROCESSED_DIR" ]] || {
                echo "[ERROR] AnswerALS RNA-seq directory not found: $ANSWERALS_PROCESSED_DIR" >&2
                return 1
            }

            [[ -s "$ANSWERALS_CHROMINFO" ]] || {
                echo "[ERROR] AnswerALS canonical chromosome sizes are missing or empty: $ANSWERALS_CHROMINFO" >&2
                return 1
            }

            while IFS= read -r donor_dir; do
                [[ -d "$donor_dir" ]] || continue
                donor="$(basename "$donor_dir")"

                # This control was excluded from the agreed cohort.
                [[ "$donor" == "CTRL-NEUEU392AE8" ]] && continue

                originals=()
                remapped_bams=()
                remapped_junctions=()

                while IFS= read -r bam; do
                    [[ -f "$bam" ]] || continue
                    case "$(basename "$bam")" in
                        *release5*|*remap*|STAR.*) continue ;;
                    esac

                    originals+=("$bam")
                    bam_name="$(basename "$bam")"
                    bam_stem="${bam_name%.sorted.bam}"
                    chimeric="$donor_dir/remap_chimeric/${bam_stem}.remap.Chimeric.out.junction"

                    if [[ -s "$chimeric" ]]; then
                        remapped_bams+=("$bam")
                        remapped_junctions+=("$chimeric")
                    fi
                done < <(find "$donor_dir" -maxdepth 1 -type f -name "*.sorted.bam" | sort)

                # Skip empty/non-sample directories.
                ((${#originals[@]} > 0)) || continue

                if ((${#remapped_bams[@]} == 0)); then
                    echo "[ERROR] No original BAM for $donor has a matching non-empty chimeric remap." >&2
                    printf '  Original BAM: %s\n' "${originals[@]}" >&2
                    return 1
                fi

                selected_index=-1
                selected_qc="not_needed"
                reason=""

                # If the donor-level link already identifies one of the
                # sample-specific remaps, preserve that selection so the
                # worker's existing donor-level circRNA outputs stay aligned.
                root_chimeric="$donor_dir/STAR.Chimeric.out.junction"
                if [[ -L "$root_chimeric" ]]; then
                    root_target="$(readlink -f "$root_chimeric" 2>/dev/null || true)"
                    [[ -n "$root_target" ]] || {
                        echo "[ERROR] Broken STAR.Chimeric.out.junction link for $donor: $root_chimeric" >&2
                        return 1
                    }

                    for i in "${!remapped_junctions[@]}"; do
                        candidate_target="$(readlink -f "${remapped_junctions[$i]}" 2>/dev/null || true)"
                        if [[ "$root_target" == "$candidate_target" ]]; then
                            selected_index="$i"
                            reason="existing STAR.Chimeric.out.junction link"
                            break
                        fi
                    done

                    if ((selected_index < 0)); then
                        echo "[ERROR] STAR.Chimeric.out.junction for $donor does not point to any eligible original BAM's remap." >&2
                        echo "  Current target: $root_target" >&2
                        printf '  Eligible remap: %s\n' "${remapped_junctions[@]}" >&2
                        return 1
                    fi
                elif [[ -e "$root_chimeric" ]]; then
                    echo "[ERROR] Expected STAR.Chimeric.out.junction to be a symlink for $donor: $root_chimeric" >&2
                    return 1
                fi

                if ((selected_index < 0)); then
                    if ((${#remapped_bams[@]} == 1)); then
                        selected_index=0
                        reason="only candidate with a matching chimeric remap"
                    else
                        # Apply the QC tiebreaker only when multiple eligible
                        # samples have remaps. Count properly paired primary
                        # alignments on canonical contigs only.
                        ensure_qc_tools

                        best_count=-1
                        for i in "${!remapped_bams[@]}"; do
                            count=$(samtools view -c \
                                -L "$QC_BED" \
                                -f 0x2 \
                                -F 0x904 \
                                "${remapped_bams[$i]}")

                            if ((count > best_count)); then
                                best_count="$count"
                                selected_index="$i"
                            fi
                        done

                        selected_qc="$best_count"
                        reason="highest canonical properly-paired primary alignment count"
                    fi
                fi

                SELECTED_BAM="${remapped_bams[$selected_index]}"
                SELECTED_CHIMERIC="${remapped_junctions[$selected_index]}"

                MASTER+=("answer_ALS|${ANSWERALS_TISSUE}|${donor}|answerals|${SELECTED_BAM}|${SELECTED_CHIMERIC}|${selected_qc}|${reason}")

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
    IFS='|' read -r EXP TISSUE SAMPLE MODE ANSWERALS_BAM ANSWERALS_CHIMERIC ANSWERALS_QC ANSWERALS_REASON <<< "$ENTRY"

    echo "=============================="
    echo "Task   : $SLURM_ARRAY_TASK_ID"
    echo "User   : $TARGET_USER"
    echo "Exp    : $EXP"
    echo "Tissue : $TISSUE"
    echo "Sample : $SAMPLE"
    echo "Mode   : $MODE"

    if [[ "$MODE" == "answerals" ]]; then
        OUTDIR="$ANSWERALS_PROCESSED_DIR/$SAMPLE"
        ROOT_CHIMERIC="$OUTDIR/STAR.Chimeric.out.junction"
        STAR_BAM="$OUTDIR/STAR.Aligned.sortedByCoord.out.bam"

        echo "BAM    : $ANSWERALS_BAM"
        echo "Remap  : $ANSWERALS_CHIMERIC"
        echo "Chosen : $ANSWERALS_REASON"
        echo "QC     : $ANSWERALS_QC"
        echo "Output : $OUTDIR"
        echo "=============================="

        [[ -s "$ANSWERALS_BAM" && -s "$ANSWERALS_CHIMERIC" ]] || {
            echo "[ERROR] Selected BAM or matching chimeric remap is missing or empty." >&2
            exit 1
        }

        # Do not let a donor-level junction link from another BAM be used.
        expected_chimeric="$(readlink -f "$ANSWERALS_CHIMERIC")"
        if [[ -L "$ROOT_CHIMERIC" ]]; then
            current_chimeric="$(readlink -f "$ROOT_CHIMERIC" 2>/dev/null || true)"
            [[ "$current_chimeric" == "$expected_chimeric" ]] || {
                echo "[ERROR] Existing STAR.Chimeric.out.junction points to a different sample for $SAMPLE." >&2
                echo "  Expected: $expected_chimeric" >&2
                echo "  Found:    $current_chimeric" >&2
                exit 1
            }
        elif [[ -e "$ROOT_CHIMERIC" ]]; then
            echo "[ERROR] Refusing to replace a non-symlink STAR.Chimeric.out.junction: $ROOT_CHIMERIC" >&2
            exit 1
        else
            ln -s "$ANSWERALS_CHIMERIC" "$ROOT_CHIMERIC"
        fi

        # If this donor already has a worker BAM link, ensure its outputs
        # correspond to the selected original BAM.
        selected_bam="$(readlink -f "$ANSWERALS_BAM")"
        if [[ -e "$STAR_BAM" || -L "$STAR_BAM" ]]; then
            existing_bam="$(readlink -f "$STAR_BAM" 2>/dev/null || true)"
            [[ "$existing_bam" == "$selected_bam" ]] || {
                echo "[ERROR] Existing STAR BAM link/output for $SAMPLE points to a different input BAM." >&2
                echo "  Selected: $selected_bam" >&2
                echo "  Existing: $existing_bam" >&2
                exit 1
            }
        fi

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
        ENTRY="${MASTER[$i]}"
        IFS='|' read -r EXP TISSUE SAMPLE MODE BAM CHIMERIC QC_COUNT REASON <<< "$ENTRY"

        if [[ "$MODE" == "answerals" ]]; then
            printf '  [%d] %s | BAM: %s | remap: %s | QC: %s | selection: %s\n' \
                "$((i + 1))" "$SAMPLE" "$BAM" "$CHIMERIC" "$QC_COUNT" "$REASON"
        else
            echo "  [$((i + 1))] $ENTRY"
        fi
    done
    echo
    echo "Would submit array: --array=1-${COUNT}%99"
    exit 0
fi

# Freeze the selected rows once so array tasks cannot independently
# recalculate or diverge from the dry-run selection.
STAMP="$(date '+%Y%m%d_%H%M%S')"
SELECTION_MANIFEST="$ANSWERALS_RNASEQ_DIR/answerALS_selected_samples_${STAMP}_$$.manifest"
printf '%s\n' "${MASTER[@]}" > "$SELECTION_MANIFEST"

echo "Saved sample selection manifest: $SELECTION_MANIFEST"
echo "Submitting array job: 1-${COUNT} (max 99 concurrent)..."
ARRAY_SUBMISSION=$(sbatch --parsable \
    --array="1-${COUNT}%99" \
    "$0" "$TARGET_USER" --manifest "$SELECTION_MANIFEST")
ARRAY_JOB_ID="${ARRAY_SUBMISSION%%;*}"
echo "Submitted array job $ARRAY_JOB_ID."
echo "Monitor with: squeue -u $TARGET_USER"
