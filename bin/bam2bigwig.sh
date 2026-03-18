#!/bin/bash
#############################################
# Author: Xianjun Dong & Zachery Wolfe
# Version: 0.9 (fixed /tmp/ collision between concurrent jobs by using job-specific temp dirs)
#############################################
set -euo pipefail

# -------- HARD-CODED PATHS (hg38) --------
CHROMINFO="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR/genome.sizes"

# use job-specific temp dir to avoid collisions between concurrent array tasks
TMPDIR=$(mktemp -d /tmp/bam2bigwig_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-$$}_XXXXXX)
trap "rm -rf $TMPDIR" EXIT

inputfile=$1   # bam or cram
split=${2:--nosplit}

[ -e "$inputfile" ] || {
    echo "Usage: bam2bigwig.sh <in.bam|cram> <-split|-nosplit>"
    exit 1
}

ext=${inputfile##*.}
bname=${inputfile%.*}

case $ext in
    bam)
        echo "Input is BAM. Converting bam → bedGraph → bigWig..."
        ;;
    cram)
        echo "Input is CRAM. Converting cram → bam → bedGraph → bigWig..."
        ;;
    *)
        echo "Unsupported input format: $ext"
        exit 1
        ;;
esac

# RPM normalization (primary reads only, unchanged logic)
RPMscale=$(bc <<< "scale=6;1000000/$(samtools view -F 0x100 -c "$inputfile")")

if [ "$split" == "-split" ]; then
    echo "bam → bw (strand-aware, fr-firststrand)"

    # PLUS transcription strand
    # plus.part1 = read2 forward: -f 0x80 -F 0x10
    # plus.part2 = read1 reverse: -f 0x40 -f 0x10
    samtools view -bh -f 0x80 -F 0x10 "$inputfile" > "$TMPDIR/plus.part1.bam"
    samtools view -bh -f 0x40 -f 0x10 "$inputfile" > "$TMPDIR/plus.part2.bam"
    samtools merge -f "$TMPDIR/plus.bam" "$TMPDIR/plus.part1.bam" "$TMPDIR/plus.part2.bam"
    bedtools genomecov -ibam "$TMPDIR/plus.bam" -bg -scale "$RPMscale" -split |
    LC_COLLATE=C sort -k1,1 -k2,2n > "$bname.plus.normalized.bedGraph"
    bedGraphToBigWig "$bname.plus.normalized.bedGraph" "$CHROMINFO" "$bname.plus.normalized.bw"

    # MINUS transcription strand
    # minus.part1 = read1 forward: -f 0x40 -F 0x10
    # minus.part2 = read2 reverse: -f 0x80 -f 0x10
    samtools view -bh -f 0x40 -F 0x10 "$inputfile" > "$TMPDIR/minus.part1.bam"
    samtools view -bh -f 0x80 -f 0x10 "$inputfile" > "$TMPDIR/minus.part2.bam"
    samtools merge -f "$TMPDIR/minus.bam" "$TMPDIR/minus.part1.bam" "$TMPDIR/minus.part2.bam"
    bedtools genomecov -ibam "$TMPDIR/minus.bam" -bg -scale "$RPMscale" -split |
    LC_COLLATE=C sort -k1,1 -k2,2n > "$bname.minus.normalized.bedGraph"
    bedGraphToBigWig "$bname.minus.normalized.bedGraph" "$CHROMINFO" "$bname.minus.normalized.bw"
fi

if [ "$split" == "-nosplit" ]; then
    echo "bam → bw (no strand split)"
    if [ "$ext" == "cram" ]; then
        samtools view -b "$inputfile" |
        bedtools genomecov -ibam stdin -bg -scale "$RPMscale" -split |
        LC_COLLATE=C sort -k1,1 -k2,2n > "$bname.normalized.bedGraph"
    else
        bedtools genomecov -ibam "$inputfile" -bg -scale "$RPMscale" -split |
        LC_COLLATE=C sort -k1,1 -k2,2n > "$bname.normalized.bedGraph"
    fi
    bedGraphToBigWig "$bname.normalized.bedGraph" "$CHROMINFO" "$bname.normalized.bw"
fi
