#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --cpus-per-task=4
#SBATCH --mem=99G
#SBATCH --time=24:00:00
#SBATCH -p day
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Submit one RNA-seq worker per sample.
# For AnswerALS, submit a cohort LeafCutter job after all sample workers succeed.

set -euo pipefail

###########################################
# Arguments
###########################################
DRY_RUN=false
LEAFCUTTER_COHORT=false
TARGET_USER=""

while (($#)); do
    case "$1" in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --leafcutter-cohort)
            LEAFCUTTER_COHORT=true
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
ANSWERALS_LEAFCUTTER_DIR="$ANSWERALS_RNASEQ_DIR/leafcutter/original_bams"

GENOME_BASE="/home/${TARGET_USER}/donglab/references/genome/Homo_sapiens/UCSC/hg38"
ANSWERALS_CHROMINFO="$GENOME_BASE/Sequence/STAR/genome.sizes.nochr"

LEAFCUTTER_BASE="/home/${TARGET_USER}/donglab/pipelines/modules/rnaseq/bin/leafcutter"

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
    HAS_ANSWERALS=false

    for entry in "${EXPERIMENTS[@]}"; do
        EXP="${entry%%|*}"
        REST="${entry#*|}"
        MODE="$REST"

        if [[ "$EXP" == "answer_ALS" && "$MODE" == "answerals" ]]; then
            HAS_ANSWERALS=true

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
# Shared AnswerALS LeafCutter cohort stage
###########################################
run_answerals_leafcutter_cohort() {
    echo "[LeafCutter] Starting AnswerALS cohort processing."

    [[ -s "$ANSWERALS_CHROMINFO" ]] || {
        echo "[ERROR] AnswerALS canonical chromosome sizes are missing or empty: $ANSWERALS_CHROMINFO" >&2
        exit 1
    }

    set +u
    source "/home/${TARGET_USER}/donglab/pipelines/modules/miniconda3/etc/profile.d/conda.sh"
    conda activate RNAseq
    set -u

    export PATH="$LEAFCUTTER_BASE/scripts:$LEAFCUTTER_BASE/clustering:$PATH"

    mkdir -p \
        "$ANSWERALS_LEAFCUTTER_DIR/juncs" \
        "$ANSWERALS_LEAFCUTTER_DIR/clusters" \
        "$ANSWERALS_LEAFCUTTER_DIR/psi"

    JUNC_LIST="$ANSWERALS_LEAFCUTTER_DIR/juncs/juncfile_list.txt"
    : > "$JUNC_LIST"

    expected=0
    for entry in "${MASTER[@]}"; do
        IFS='|' read -r EXP TISSUE SAMPLE MODE BAM <<< "$entry"
        [[ "$MODE" == "answerals" ]] || continue
        ((expected += 1))

        OUTDIR="$ANSWERALS_PROCESSED_DIR/$SAMPLE"
        JUNC="$OUTDIR/leafcutter/juncs/${SAMPLE}.junc"

        [[ -f "$OUTDIR/.status.RNAseq.leafcutter_junctions" && -s "$JUNC" ]] || {
            echo "[ERROR] Missing or incomplete AnswerALS junction output for $SAMPLE: $JUNC" >&2
            exit 1
        }

        # Enforce the canonical-contig rule again at cohort input. This also
        # cleans up junction files left by earlier runs before the worker
        # filtered them. Replace each derived .junc file atomically.
        FILTERED_JUNC="$(mktemp "${JUNC}.tmp.XXXXXX")"
        awk '
            NR == FNR { allowed[$1] = 1; next }
            NF && ($1 in allowed)
        ' "$ANSWERALS_CHROMINFO" "$JUNC" > "$FILTERED_JUNC"

        [[ -s "$FILTERED_JUNC" ]] || {
            rm -f "$FILTERED_JUNC"
            echo "[ERROR] No canonical-contig LeafCutter junctions remain for $SAMPLE: $JUNC" >&2
            exit 1
        }

        mv "$FILTERED_JUNC" "$JUNC"
        printf '%s\n' "$JUNC" >> "$JUNC_LIST"
    done

    actual=$(wc -l < "$JUNC_LIST")
    if [[ "$actual" -ne "$expected" || "$actual" -eq 0 ]]; then
        echo "[ERROR] Expected $expected junction files, found $actual." >&2
        exit 1
    fi

    CLUSTER_DIR="$ANSWERALS_LEAFCUTTER_DIR/clusters"
    PSI_DIR="$ANSWERALS_LEAFCUTTER_DIR/psi"

    # Create the cohort's refined clusters from all per-sample junction files.
    python "$LEAFCUTTER_BASE/clustering/leafcutter_cluster.py" \
        -j "$JUNC_LIST" \
        -o leafcutter_clusters \
        -r "$CLUSTER_DIR"

    # Quantify using the refined clusters just created above.
    python "$LEAFCUTTER_BASE/clustering/leafcutter_quant_only.py" \
        -j "$JUNC_LIST" \
        -o leafcutter_clusters \
        -r "$CLUSTER_DIR"

    CLUSTERS="$CLUSTER_DIR/leafcutter_clusters_refined"
    [[ -s "$CLUSTERS" ]] || {
        echo "[ERROR] LeafCutter refined cluster file was not created: $CLUSTERS" >&2
        exit 1
    }

    # Produce one PSI file per sample from the shared cohort clusters.
    for entry in "${MASTER[@]}"; do
        IFS='|' read -r EXP TISSUE SAMPLE MODE BAM <<< "$entry"
        [[ "$MODE" == "answerals" ]] || continue

        JUNC="$ANSWERALS_PROCESSED_DIR/$SAMPLE/leafcutter/juncs/${SAMPLE}.junc"
        PSI="$PSI_DIR/${SAMPLE}.leafcutter.PSI.tsv"

        python "$LEAFCUTTER_BASE/scripts/leafcutter_psi.py" \
            "$CLUSTERS" "$JUNC" "$PSI"

        {
            head -n 1 "$PSI"
            tail -n +2 "$PSI" | sort -k7,7nr
        } > "$PSI_DIR/${SAMPLE}.leafcutter.PSI.sorted.tsv"
    done

    touch "$ANSWERALS_LEAFCUTTER_DIR/.status.leafcutter_cohort"
    echo "[LeafCutter] Cohort processing completed for $actual samples."
}

###########################################
# Cohort-only job mode
###########################################
if "$LEAFCUTTER_COHORT"; then
    build_master
    run_answerals_leafcutter_cohort
    exit 0
fi

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

    if "$HAS_ANSWERALS"; then
        echo "Would submit AnswerALS LeafCutter cohort job after the array succeeds."
    fi
    exit 0
fi

echo "Submitting array job: 1-${COUNT} (max 99 concurrent)..."
ARRAY_SUBMISSION=$(sbatch --parsable --array="1-${COUNT}%99" "$0" "$TARGET_USER")
ARRAY_JOB_ID="${ARRAY_SUBMISSION%%;*}"
echo "Submitted array job $ARRAY_JOB_ID."

if "$HAS_ANSWERALS"; then
    echo "Submitting AnswerALS LeafCutter cohort job, dependent on successful array completion..."
    COHORT_JOB_ID=$(sbatch --parsable \
        --dependency="afterok:${ARRAY_JOB_ID}" \
        "$0" --leafcutter-cohort "$TARGET_USER")
    echo "Submitted LeafCutter cohort job ${COHORT_JOB_ID%%;*}."
fi

echo "Monitor with: squeue -u $TARGET_USER"
