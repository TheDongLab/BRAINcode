#!/bin/bash
#############################################
# Author: Xianjun Dong & Zachery Wolfe
# Version: 0.2 (patched for Zach's updated RNAseq pipeline)
#############################################

set -euo pipefail

# ----------------------------
# Inputs
# ----------------------------
inputfile="$1"                     # BAM/SAM/BED
split="${2:-"-split"}"             # "-split" or "-nosplit"
normalizationFactor="${3:-0}"      # if 0, no normalization
region="${4:-}"                    # optional: "chrX", "chrX:1000-2000", etc.

# ----------------------------
# Paths (match RNAseq pipeline)
# ----------------------------
GENOME_DIR=/home/zw529/donglab/pipelines/genome
ANNOTATION="$GENOME_DIR"                     # ChromInfo.txt lives here
SAMPLE_DIR="$(dirname "$inputfile")"
TMPDIR="$SAMPLE_DIR/tmp"
mkdir -p "$TMPDIR"
export TMPDIR

bname=$(basename "$inputfile" .bam)
bname="${bname%.sam}"
bname="${bname%.bed}"

species=hg38   # match your genome

# ----------------------------
# Check input
# ----------------------------
if [ ! -e "$inputfile" ]; then
    echo "ERROR: input file $inputfile not found!"
    exit 1
fi

# ----------------------------
# Convert input to BED12
# ----------------------------
ext="${inputfile##*.}"

case "$ext" in
    bam|BAM|Bam)
        echo "Input is BAM. Converting bam -> bed..."
        samtools view "$inputfile" $region | sam2bed -v bed12=T -v sCol=NH > "$SAMPLE_DIR/$bname.bed";;
    sam|SAM|Sam)
        echo "Input is SAM. Converting sam -> bed..."
        sam2bed -v bed12=T -v sCol=NH "$inputfile" > "$SAMPLE_DIR/$bname.bed";;
    bed|BED|Bed|bed6|bed12)
        echo "Input is BED. Using as-is."
        cp "$inputfile" "$SAMPLE_DIR/$bname.bed";;
    *)
        echo "Unsupported input format: $ext"
        exit 1;;
esac

BED="$SAMPLE_DIR/$bname.bed"

# ----------------------------
# Split by strand if requested
# ----------------------------
if [ "$split" == "-split" ]; then
    echo "Splitting BED by strand..."
    > "$SAMPLE_DIR/$bname.bed+"
    > "$SAMPLE_DIR/$bname.bed-"
    awk -v OUT="$SAMPLE_DIR/$bname" '{print > OUT".bed"$6}' "$BED"

    ncol=$(awk 'END{print NF}' "$BED")
    if [ "$ncol" -eq 12 ]; then
        # bed12
        sort -k1,1 "$SAMPLE_DIR/$bname.bed+" | \
        bedItemOverlapCount "$species" -bed12 -chromSize="$ANNOTATION/ChromInfo.txt" stdin | \
        sort -k1,1 -k2,2n > "$SAMPLE_DIR/$bname.plus.bedGraph"

        sort -k1,1 "$SAMPLE_DIR/$bname.bed-" | \
        bedItemOverlapCount "$species" -bed12 -chromSize="$ANNOTATION/ChromInfo.txt" stdin | \
        sort -k1,1 -k2,2n | \
        awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > "$SAMPLE_DIR/$bname.minus.bedGraph"
    else
        # bed6
        sort -k1,1 "$SAMPLE_DIR/$bname.bed+" | \
        bedItemOverlapCount "$species" -chromSize="$ANNOTATION/ChromInfo.txt" stdin | \
        sort -k1,1 -k2,2n > "$SAMPLE_DIR/$bname.plus.bedGraph"

        sort -k1,1 "$SAMPLE_DIR/$bname.bed-" | \
        bedItemOverlapCount "$species" -chromSize="$ANNOTATION/ChromInfo.txt" stdin | \
        sort -k1,1 -k2,2n | \
        awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > "$SAMPLE_DIR/$bname.minus.bedGraph"
    fi

    # Convert to BigWig
    if [ "$normalizationFactor" -eq 0 ]; then
        bedGraphToBigWig "$SAMPLE_DIR/$bname.plus.bedGraph" "$ANNOTATION/ChromInfo.txt" "$SAMPLE_DIR/$bname.plus.bw"
        bedGraphToBigWig "$SAMPLE_DIR/$bname.minus.bedGraph" "$ANNOTATION/ChromInfo.txt" "$SAMPLE_DIR/$bname.minus.bw"
    else
        awk -v tmr="$normalizationFactor" 'BEGIN{OFS="\t"}{$4=$4*1e6/tmr; print}' "$SAMPLE_DIR/$bname.plus.bedGraph" > "$SAMPLE_DIR/$bname.plus.normalized.bedGraph"
        awk -v tmr="$normalizationFactor" 'BEGIN{OFS="\t"}{$4=$4*1e6/tmr; print}' "$SAMPLE_DIR/$bname.minus.bedGraph" > "$SAMPLE_DIR/$bname.minus.normalized.bedGraph"
        bedGraphToBigWig "$SAMPLE_DIR/$bname.plus.normalized.bedGraph" "$ANNOTATION/ChromInfo.txt" "$SAMPLE_DIR/$bname.plus.normalized.bw"
        bedGraphToBigWig "$SAMPLE_DIR/$bname.minus.normalized.bedGraph" "$ANNOTATION/ChromInfo.txt" "$SAMPLE_DIR/$bname.minus.normalized.bw"
    fi
fi

# ----------------------------
# No split
# ----------------------------
if [ "$split" == "-nosplit" ]; then
    echo "Processing BED without strand split..."
    ncol=$(awk 'END{print NF}' "$BED")
    if [ "$ncol" -eq 12 ]; then
        sort -k1,1 "$BED" | \
        bedItemOverlapCount "$species" -bed12 -chromSize="$ANNOTATION/ChromInfo.txt" stdin | \
        sort -k1,1 -k2,2n > "$SAMPLE_DIR/$bname.bedGraph"
    else
        sort -k1,1 "$BED" | \
        bedItemOverlapCount "$species" -chromSize="$ANNOTATION/ChromInfo.txt" stdin | \
        sort -k1,1 -k2,2n > "$SAMPLE_DIR/$bname.bedGraph"
    fi

    if [ "$normalizationFactor" -eq 0 ]; then
        bedGraphToBigWig "$SAMPLE_DIR/$bname.bedGraph" "$ANNOTATION/ChromInfo.txt" "$SAMPLE_DIR/$bname.bw"
    else
        awk -v tmr="$normalizationFactor" 'BEGIN{OFS="\t"}{$4=$4*1e6/tmr; print}' "$SAMPLE_DIR/$bname.bedGraph" > "$SAMPLE_DIR/$bname.normalized.bedGraph"
        bedGraphToBigWig "$SAMPLE_DIR/$bname.normalized.bedGraph" "$ANNOTATION/ChromInfo.txt" "$SAMPLE_DIR/$bname.normalized.bw"
    fi
fi

echo "bam2bigwig finished: $SAMPLE_DIR/$bname.*.bw"
