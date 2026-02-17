#!/bin/bash
###########################################
# _bam2annotation.sh
# Use safe annotation folder
###########################################

set -euo pipefail

BAM_FILE=$1
OUT_PREFIX=$2
GENOME_DIR=/home/zw529/donglab/pipelines/genome/annotation_safe
CHROMSIZES="$GENOME_DIR/chromsizes.safe.txt"

BED_DIR="${OUT_PREFIX}_beds"
mkdir -p "$BED_DIR"

BAM_BED="$BED_DIR/$(basename "$BAM_FILE" .bam).bed"
bedtools bamtobed -bed12 -i "$BAM_FILE" > "$BAM_BED"

# Intersections with safe BEDs
for f in exons introns genes genic intergenic intergenic_notneargene rRNA non_rRNA_mt; do
    bedtools intersect -s -wa -a "$BAM_BED" -b "$GENOME_DIR/$f.sorted.bed" > "$BED_DIR/$f.bed"
done

# mtRNA
grep -P '^chrM\t' "$BAM_BED" > "$BED_DIR/mtRNA.bed"

# Summary
SUMMARY="${OUT_PREFIX}.summary.txt"
{
  echo -e "category\tcount"
  for f in exons introns rRNA mtRNA LINE SINE intergenic intergenic_notneargene; do
      case $f in
          LINE|SINE) echo -e "$f\t0" ;;
          *) echo -e "$f\t$(wc -l < "$BED_DIR/$f.bed")" ;;
      esac
  done
} > "$SUMMARY"

echo "[`date`] Annotation summary written to $SUMMARY"
