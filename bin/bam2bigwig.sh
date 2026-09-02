#!/bin/bash
#############################################
# Author: Xianjun Dong & Zachery Wolfe
# Version: 1.0
# Strand-specific BigWig splitting for targetALS reverse-stranded / fr-firststrand libraries.
#############################################
set -euo pipefail

# ============================================================
# HARD-CODED PATHS (hg38)
# ============================================================

CHROMINFO="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR/genome.sizes"

# ============================================================
# TEMP DIRECTORY
# ============================================================

TMPDIR=$(mktemp -d /tmp/bam2bigwig_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-$$}_XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

# ============================================================
# INPUTS
# ============================================================

inputfile=$1
split=${2:--nosplit}

if [ ! -e "$inputfile" ]; then
    echo "Usage: bam2bigwig.sh <in.bam|cram> <-split|-nosplit>"
    exit 1
fi

ext=${inputfile##*.}
bname=${inputfile%.*}

case "$ext" in
    bam)
        echo "Input is BAM. Converting BAM → bedGraph → BigWig..."
        ;;
    cram)
        echo "Input is CRAM. Converting CRAM → bedGraph → BigWig..."
        ;;
    *)
        echo "ERROR: Unsupported input format: $ext"
        exit 1
        ;;
esac

# ============================================================
# RPM NORMALIZATION
#
# Keep original normalization logic:
#   1e6 / number of primary alignments
#
# Remove secondary alignments:
#   0x100
# ============================================================

PRIMARY_READS=$(samtools view -F 0x100 -c "$inputfile")

if [ "$PRIMARY_READS" -eq 0 ]; then
    echo "ERROR: No primary alignments found in $inputfile"
    exit 1
fi

RPMscale=$(bc <<< "scale=6;1000000/$PRIMARY_READS")

echo "Primary alignments : $PRIMARY_READS"
echo "RPM scale          : $RPMscale"

# ============================================================
# STRAND-SPECIFIC BIGWIG
# ============================================================

if [ "$split" == "-split" ]; then

    echo "bam → bw (strand-aware, fr-firststrand / reverse-stranded)"
    echo
    echo "Strand definitions:"
    echo "  PLUS  = R2 forward + R1 reverse"
    echo "  MINUS = R1 forward + R2 reverse"
    echo

    # ========================================================
    # PLUS TRANSCRIPTION STRAND
    #
    # targetALS is reverse-stranded / fr-firststrand:
    #
    #   + transcript:
    #       R1 = reverse
    #       R2 = forward
    #
    # PLUS = R2 forward + R1 reverse
    # ========================================================

    echo "[PLUS] Extracting R2 forward..."
    samtools view \
        -bh \
        -f 0x80 \
        -F 0x10 \
        "$inputfile" \
        > "$TMPDIR/plus.part1.bam"

    echo "[PLUS] Extracting R1 reverse..."
    samtools view \
        -bh \
        -f 0x40 \
        -f 0x10 \
        "$inputfile" \
        > "$TMPDIR/plus.part2.bam"

    echo "[PLUS] Merging orientation classes..."
    samtools merge \
        -f \
        "$TMPDIR/plus.bam" \
        "$TMPDIR/plus.part1.bam" \
        "$TMPDIR/plus.part2.bam"

    echo "[PLUS] Generating normalized bedGraph..."
    bedtools genomecov \
        -ibam "$TMPDIR/plus.bam" \
        -bg \
        -scale "$RPMscale" \
        -split |
    LC_COLLATE=C sort -k1,1 -k2,2n \
        > "$bname.plus.normalized.bedGraph"

    echo "[PLUS] Generating BigWig..."
    bedGraphToBigWig \
        "$bname.plus.normalized.bedGraph" \
        "$CHROMINFO" \
        "$bname.plus.normalized.bw"

    # ========================================================
    # MINUS TRANSCRIPTION STRAND
    #
    # targetALS is reverse-stranded / fr-firststrand:
    #
    #   - transcript:
    #       R1 = forward
    #       R2 = reverse
    #
    # MINUS = R1 forward + R2 reverse
    # ========================================================

    echo "[MINUS] Extracting R1 forward..."
    samtools view \
        -bh \
        -f 0x40 \
        -F 0x10 \
        "$inputfile" \
        > "$TMPDIR/minus.part1.bam"

    echo "[MINUS] Extracting R2 reverse..."
    samtools view \
        -bh \
        -f 0x80 \
        -f 0x10 \
        "$inputfile" \
        > "$TMPDIR/minus.part2.bam"

    echo "[MINUS] Merging orientation classes..."
    samtools merge \
        -f \
        "$TMPDIR/minus.bam" \
        "$TMPDIR/minus.part1.bam" \
        "$TMPDIR/minus.part2.bam"

    echo "[MINUS] Generating normalized bedGraph..."
    bedtools genomecov \
        -ibam "$TMPDIR/minus.bam" \
        -bg \
        -scale "$RPMscale" \
        -split |
    LC_COLLATE=C sort -k1,1 -k2,2n \
        > "$bname.minus.normalized.bedGraph"

    echo "[MINUS] Generating BigWig..."
    bedGraphToBigWig \
        "$bname.minus.normalized.bedGraph" \
        "$CHROMINFO" \
        "$bname.minus.normalized.bw"

    echo
    echo "Strand-specific BigWigs completed:"
    echo "  $bname.plus.normalized.bw"
    echo "  $bname.minus.normalized.bw"

# ============================================================
# NON-STRAND-SPECIFIC BIGWIG
# ============================================================

elif [ "$split" == "-nosplit" ]; then

    echo "bam → bw (no strand split)"

    if [ "$ext" == "cram" ]; then

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

    echo
    echo "BigWig completed:"
    echo "  $bname.normalized.bw"

else

    echo "ERROR: Second argument must be -split or -nosplit"
    exit 1

fi
