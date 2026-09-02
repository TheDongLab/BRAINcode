#!/bin/bash
#SBATCH --job-name=infer_strand
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/infer_strand.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/infer_strand.err

set -euo pipefail

set +u
source /home/zw529/donglab/pipelines/modules/miniconda3/etc/profile.d/conda.sh
conda activate RNAseq
set -u

command -v infer_experiment.py >/dev/null || {
    echo "ERROR: infer_experiment.py not found in RNAseq environment"
    exit 1
}

ROOT="/home/zw529/donglab/data/target_ALS"
BED="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.labeled.transcript.bed12"
OUT="$ROOT/QTL/infer_experiment_all"
LIST="$OUT/samples.tsv"
mkdir -p "$OUT/results"

# ------------------------------------------------------------
# INITIAL CALL: build sample list + submit array
# ------------------------------------------------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo -e "index\ttissue\tsample\tdir" > "$LIST"
    i=0

    while IFS= read -r D; do
        [[ -d "$D" ]] || continue
        ((++i))
        tissue=$(basename "$(dirname "$(dirname "$(dirname "$D")")")")
        sample=$(basename "$D")
        echo -e "$i\t$tissue\t$sample\t$D" >> "$LIST"
    done < <(find "$ROOT" -type d -path '*/RNAseq/Processed' -print0 | while IFS= read -r -d '' P; do find "$P" -mindepth 1 -maxdepth 1 -type d; done | sort)

    echo "Samples found: $i"
    sbatch --array="1-${i}%500" "$0"
    exit 0
fi

# ------------------------------------------------------------
# ARRAY TASK
# ------------------------------------------------------------
LINE=$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "$LIST")
IFS=$'\t' read -r IDX TISSUE SAMPLE D <<< "$LINE"

RESULT="$OUT/results/${IDX}.tsv"
LOG="$OUT/results/${IDX}.infer.txt"

echo "[$IDX] $TISSUE | $SAMPLE"

# Prefer the genuine targetALS source BAM, NOT circ/remap test BAMs.
BAM=$(find "$D" -maxdepth 1 -type f -name '*.final.bam' \
      ! -name '*remap*' ! -name '*circ_boundary*' | head -1 || true)

# Fallback to STAR BAM only if its resolved target is sensible.
if [[ -z "$BAM" && -e "$D/STAR.Aligned.sortedByCoord.out.bam" ]]; then
    CAND=$(readlink -f "$D/STAR.Aligned.sortedByCoord.out.bam" || true)
    if [[ -n "$CAND" && -f "$CAND" && "$CAND" != *remap* && "$CAND" != *circ_boundary* ]]; then
        BAM="$CAND"
    fi
fi

if [[ -z "$BAM" ]]; then
    echo -e "$IDX\t$TISSUE\t$SAMPLE\tNO_VALID_BAM\tNA\tNA\tNA\tNO_BAM" > "$RESULT"
    exit 0
fi

infer_experiment.py -i "$BAM" -r "$BED" -s 200000 > "$LOG" 2>&1 || {
    echo -e "$IDX\t$TISSUE\t$SAMPLE\t$BAM\tNA\tNA\tNA\tFAILED" > "$RESULT"
    exit 0
}

FAILED=$(awk -F': ' '/Fraction of reads failed to determine/{print $2}' "$LOG")
P1=$(awk -F': ' '/1\+\+,1--,2\+-,2-\+/{print $2}' "$LOG")
P2=$(awk -F': ' '/1\+-,1-\+,2\+\+,2--/{print $2}' "$LOG")

STATUS=$(awk -v p1="$P1" -v p2="$P2" 'BEGIN{
    if(p1=="" || p2=="") print "PARSE_FAILED";
    else if(p2 >= 5*p1) print "FIRSTSTRAND";
    else if(p1 >= 5*p2) print "SECONDSTRAND";
    else print "AMBIGUOUS";
}')

echo -e "$IDX\t$TISSUE\t$SAMPLE\t$BAM\t$FAILED\t$P1\t$P2\t$STATUS" > "$RESULT"
