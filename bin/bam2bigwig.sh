#!/bin/bash
#############################################
# Strand-specific BigWig generation.
# Default reference supports existing TargetALS calls.
# An optional third argument supplies AnswerALS chromosome sizes.
#############################################

set -euo pipefail

DEFAULT_CHROMINFO="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR/genome.sizes"

###########################################
# Inputs
###########################################
inputfile="${1:-}"
split="${2:--nosplit}"
CHROMINFO="${3:-$DEFAULT_CHROMINFO}"

if [[ -z "$inputfile" || ! -e "$inputfile" ]]; then
    echo "Usage: bam2bigwig.sh <in.bam|cram> <-split|-nosplit> [genome.sizes]" >&2
    exit 1
fi

if [[ ! -s "$CHROMINFO" ]]; then
    echo "[ERROR] Chromosome sizes file is missing or empty: $CHROMINFO" >&2
    exit 1
fi

ext="${inputfile##*.}"
bname="${inputfile%.*}"

case "$ext" in
    bam|cram)
        echo "Input is ${ext^^}. Generating normalized BigWig..."
        ;;
    *)
        echo "[ERROR] Unsupported input format: $ext" >&2
        exit 1
        ;;
esac

case "$split" in
    -split|-nosplit)
        ;;
    *)
        echo "[ERROR] Second argument must be -split or -nosplit." >&2
        exit 1
        ;;
esac

###########################################
# Check that BAM/CRAM references occur in the supplied sizes file.
# This catches a chr-prefixed/unprefixed mismatch before conversion.
###########################################
while IFS=$'\t' read -r chrom length rest; do
    [[ -n "$chrom" && -n "$length" ]] || continue
    if ! awk -F'\t' -v chrom="$chrom" '$1 == chrom { found=1; exit } END { exit !found }' "$CHROMINFO"; then
        echo "[ERROR] BAM/CRAM contig '$chrom' is absent from $CHROMINFO" >&2
        exit 1
    fi
done < <(
    samtools view -H "$inputfile" |
        awk -F'\t' '$1 == "@SQ" {
            sn=""; ln="";
            for (i=2; i<=NF; i++) {
                if ($i ~ /^SN:/) sn=substr($i,4);
                if ($i ~ /^LN:/) ln=substr($i,4);
            }
            if (sn != "" && ln != "") print sn "\t" ln;
        }'
)

###########################################
# Temporary directory
###########################################
WORK_TMPDIR=$(mktemp -d "/tmp/bam2bigwig_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-$$}_XXXXXX")
trap 'rm -rf "$WORK_TMPDIR"' EXIT

###########################################
# RPM normalization: primary alignments only
###########################################
PRIMARY_READS=$(samtools view -F 0x100 -c "$inputfile")

if [[ "$PRIMARY_READS" -eq 0 ]]; then
    echo "[ERROR] No primary alignments found in $inputfile" >&2
    exit 1
fi

RPMscale=$(bc <<< "scale=6;1000000/$PRIMARY_READS")

echo "Primary alignments : $PRIMARY_READS"
echo "RPM scale          : $RPMscale"
echo "Chromosome sizes   : $CHROMINFO"

###########################################
# Strand-specific BigWigs
#
# For reverse-stranded / fr-firststrand libraries:
#   PLUS  = R2 forward + R1 reverse
#   MINUS = R1 forward + R2 reverse
###########################################
if [[ "$split" == "-split" ]]; then
    echo "Generating strand-specific BigWigs (reverse-stranded / fr-firststrand)."

    echo "[PLUS] Extracting R2 forward reads..."
    samtools view -bh -f 0x80 -F 0x10 "$inputfile" \
        > "$WORK_TMPDIR/plus.part1.bam"

    echo "[PLUS] Extracting R1 reverse reads..."
    samtools view -bh -f 0x40 -f 0x10 "$inputfile" \
        > "$WORK_TMPDIR/plus.part2.bam"

    echo "[PLUS] Merging read orientations..."
    samtools merge -f "$WORK_TMPDIR/plus.bam" \
        "$WORK_TMPDIR/plus.part1.bam" \
        "$WORK_TMPDIR/plus.part2.bam"

    echo "[PLUS] Creating bedGraph..."
    bedtools genomecov \
        -ibam "$WORK_TMPDIR/plus.bam" \
        -bg \
        -scale "$RPMscale" \
        -split |
        LC_COLLATE=C sort -k1,1 -k2,2n \
        > "$bname.plus.normalized.bedGraph"

    echo "[PLUS] Creating BigWig..."
    bedGraphToBigWig \
        "$bname.plus.normalized.bedGraph" \
        "$CHROMINFO" \
        "$bname.plus.normalized.bw"

    echo "[MINUS] Extracting R1 forward reads..."
    samtools view -bh -f 0x40 -F 0x10 "$inputfile" \
        > "$WORK_TMPDIR/minus.part1.bam"

    echo "[MINUS] Extracting R2 reverse reads..."
    samtools view -bh -f 0x80 -f 0x10 "$inputfile" \
        > "$WORK_TMPDIR/minus.part2.bam"

    echo "[MINUS] Merging read orientations..."
    samtools merge -f "$WORK_TMPDIR/minus.bam" \
        "$WORK_TMPDIR/minus.part1.bam" \
        "$WORK_TMPDIR/minus.part2.bam"

    echo "[MINUS] Creating bedGraph..."
    bedtools genomecov \
        -ibam "$WORK_TMPDIR/minus.bam" \
        -bg \
        -scale "$RPMscale" \
        -split |
        LC_COLLATE=C sort -k1,1 -k2,2n \
        > "$bname.minus.normalized.bedGraph"

    echo "[MINUS] Creating BigWig..."
    bedGraphToBigWig \
        "$bname.minus.normalized.bedGraph" \
        "$CHROMINFO" \
        "$bname.minus.normalized.bw"

    echo "Strand-specific BigWigs completed:"
    echo "  $bname.plus.normalized.bw"
    echo "  $bname.minus.normalized.bw"

###########################################
# Unsplit BigWig
###########################################
else
    echo "Generating unsplit BigWig."

    if [[ "$ext" == "cram" ]]; then
        samtools view -b "$inputfile" |
            bedtools genomecov \
                -ibam stdin \
                -bg \
                -scale "$RPMscale" \
                -split |
            LC_COLLATE=C sort -k1,1 -k2,2n \
            > "$bname.normalized.bedGraph"
    else
        bedtools genomecov \
            -ibam "$inputfile" \
            -bg \
            -scale "$RPMscale" \
            -split |
            LC_COLLATE=C sort -k1,1 -k2,2n \
            > "$bname.normalized.bedGraph"
    fi

    bedGraphToBigWig \
        "$bname.normalized.bedGraph" \
        "$CHROMINFO" \
        "$bname.normalized.bw"

    echo "BigWig completed:"
    echo "  $bname.normalized.bw"
fi
