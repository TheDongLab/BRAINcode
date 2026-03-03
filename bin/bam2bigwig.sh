#!/bin/bash
#############################################
# Author: Xianjun Dong & Zachery Wolfe
# Email: xdong@rics.bwh.harvard.edu 
# Version: 0.8

# Fixed bam2bigwig.sh
# - removes neurogen TMPDIR
# - removes config file dependency
# - hardcodes hg38 ChromInfo path
# - keeps original logic intact
# - adds strand logic (i.e., Gene "+"  =  (read2 & forward)  OR  (read1 & reverse) and Gene "-"  =  (read1 & forward)  OR  (read2 & reverse)) for + and - bigWig creation
#############################################

set -euo pipefail

# -------- HARD-CODED PATHS (hg38) --------
CHROMINFO="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR/genome.sizes"

# use system temp instead of non-existent /data/neurogen
export TMPDIR=/tmp

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
    samtools view -bh \
        -f 0x80 -F 0x10 "$inputfile" \
        "$inputfile" |
    samtools view -bh \
        -f 0x40 -f 0x10 "$inputfile" > /tmp/plus.part2.bam

    samtools view -bh -f 0x80 -F 0x10 "$inputfile" > /tmp/plus.part1.bam
    samtools view -bh -f 0x40 -f 0x10 "$inputfile" > /tmp/plus.part2.bam

    samtools merge -f /tmp/plus.bam /tmp/plus.part1.bam /tmp/plus.part2.bam

    bedtools genomecov -ibam /tmp/plus.bam -bg -scale "$RPMscale" -split |
    LC_COLLATE=C sort -k1,1 -k2,2n > "$bname.plus.normalized.bedGraph"

    bedGraphToBigWig "$bname.plus.normalized.bedGraph" "$CHROMINFO" "$bname.plus.normalized.bw"

    # MINUS transcription strand
    samtools view -bh -f 0x40 -F 0x10 "$inputfile" > /tmp/minus.part1.bam
    samtools view -bh -f 0x80 -f 0x10 "$inputfile" > /tmp/minus.part2.bam

    samtools merge -f /tmp/minus.bam /tmp/minus.part1.bam /tmp/minus.part2.bam

    bedtools genomecov -ibam /tmp/minus.bam -bg -scale "$RPMscale" -split |
    LC_COLLATE=C sort -k1,1 -k2,2n > "$bname.minus.normalized.bedGraph"

    bedGraphToBigWig "$bname.minus.normalized.bedGraph" "$CHROMINFO" "$bname.minus.normalized.bw"

    rm -f /tmp/plus.* /tmp/minus.*
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
