#!/bin/bash
#############################################
# Author: Xianjun Dong & Zachery Wolfe
# Version: 0.7 (deeptools-based, strand-aware, auto-indexing)
# bam2bigwig — robust, normalized, chr-safe, tmp cleanup
#############################################

set -euo pipefail

# ----------------------------
# Genome / chrom sizes (hg38 STAR)
# ----------------------------
GENOME_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR
CHROMSIZES="$GENOME_DIR/genome.sizes"

# ----------------------------
# Hard requirements
# ----------------------------
req=(samtools bedtools awk sort module)
for x in "${req[@]}"; do
    command -v "$x" >/dev/null || { echo "MISSING: $x"; exit 1; }
done

module load deepTools

# ----------------------------
# Inputs
# ----------------------------
inputfile="$1"         # BAM, SAM, or BED
split="${2:--split}"   # -split or -nosplit
norm_reads="${3:-0}"   # optional: primary reads for normalization
region=""              # optional region can remain empty

SAMPLE_DIR="$(cd "$(dirname "$inputfile")" && pwd)"
TMPDIR="$SAMPLE_DIR/tmp"
mkdir -p "$TMPDIR"

bname="$(basename "$inputfile")"
bname="${bname%.bam}"
bname="${bname%.sam}"
bname="${bname%.bed}"

export LC_COLLATE=C
export TMPDIR

# ----------------------------
# Convert input -> BED12 if needed
# ----------------------------
BED="$SAMPLE_DIR/$bname.bed"

if [[ "${inputfile##*.}" =~ ^bam|BAM$ ]]; then
    samtools view -h "$inputfile" | bedtools bamtobed -bed12 -i stdin > "$BED"
elif [[ "${inputfile##*.}" =~ ^sam|SAM$ ]]; then
    samtools view -h "$inputfile" | bedtools bamtobed -bed12 -i stdin > "$BED"
elif [[ "${inputfile##*.}" =~ ^bed|BED$ ]]; then
    [[ "$inputfile" != "$BED" ]] && cp "$inputfile" "$BED"
else
    echo "Unsupported input: $inputfile"
    exit 1
fi

[[ -s "$BED" ]] || { echo "BED empty after conversion"; exit 1; }

# ----------------------------
# Normalize chr naming
# ----------------------------
chr_mode=$(awk 'NR==1{print ($1 ~ /^chr/)}' "$CHROMSIZES")
if [[ "$chr_mode" == 1 ]]; then
    awk 'BEGIN{OFS="\t"}{$1=($1 ~ /^chr/)?$1:"chr"$1; print}' "$BED" > "$BED.tmp"
else
    awk 'BEGIN{OFS="\t"}{$1=gensub(/^chr/,"","",$1); print}' "$BED" > "$BED.tmp"
fi
mv "$BED.tmp" "$BED"

# ----------------------------
# Sort BED
# ----------------------------
sort -k1,1 -k2,2n "$BED" > "$BED.sorted"
mv "$BED.sorted" "$BED"

# ----------------------------
# Coverage helper via deepTools
# ----------------------------
bam_to_bw () {
    local infile="$1"
    local outfile="$2"

    [[ -f "$infile.bai" ]] || samtools index "$infile"

    bamCoverage -b "$infile" -o "$outfile" \
        --binSize 50 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --ignoreDuplicates \
        --extendReads 150 \
        --numberOfProcessors 8
}

# ----------------------------
# Split by strand if requested
# ----------------------------
if [[ "$split" == "-split" ]]; then
    awk '$6=="+"' "$BED" > "$BED.plus"
    awk '$6=="-"' "$BED" > "$BED.minus"

    BEDTOOLS_BAM_PLUS="$TMPDIR/$bname.plus.bam"
    BEDTOOLS_BAM_MINUS="$TMPDIR/$bname.minus.bam"
    bedtools bedtobam -i "$BED.plus" -g "$CHROMSIZES" > "$BEDTOOLS_BAM_PLUS"
    bedtools bedtobam -i "$BED.minus" -g "$CHROMSIZES" > "$BEDTOOLS_BAM_MINUS"

    bam_to_bw "$BEDTOOLS_BAM_PLUS" "$SAMPLE_DIR/$bname.plus.bw"
    bam_to_bw "$BEDTOOLS_BAM_MINUS" "$SAMPLE_DIR/$bname.minus.bw"
else
    # Use original BAM if available, else convert BED -> BAM
    if [[ "${inputfile##*.}" =~ ^bam|BAM$ ]]; then
        bam_to_bw "$inputfile" "$SAMPLE_DIR/$bname.bw"
    else
        BEDTOOLS_BAM="$TMPDIR/$bname.tmp.bam"
        bedtools bedtobam -i "$BED" -g "$CHROMSIZES" > "$BEDTOOLS_BAM"
        bam_to_bw "$BEDTOOLS_BAM" "$SAMPLE_DIR/$bname.bw"
    fi
fi

rm -rf "$TMPDIR"
echo "OK: $bname bigWig(s) created"
