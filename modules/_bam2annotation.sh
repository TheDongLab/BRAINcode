#!/bin/bash
# _bam2annotation.sh: simplified but fully compatible read-counts for _bam2annotation.R
# usage: ./_bam2annotation.sh accepted_hits.bam

if [ $# -ne 1 ]; then
    echo "Usage: $0 accepted_hits.bam"
    exit 1
fi

INPUT_BAM=$1
ANNOTATION=/home/zw529/donglab/pipelines/genome
GENOME_BED12=$ANNOTATION/hg38.bed12   # generated from GTF
RRNA_BED=$ANNOTATION/rRNA.bed
MT_BED=$ANNOTATION/chrM.bed

test -e $INPUT_BAM || { echo "$INPUT_BAM not found"; exit 1; }

# convert BAM to BED once
BAM_BED=$INPUT_BAM.bed
bedtools bamtobed -i $INPUT_BAM > $BAM_BED

# total reads
TOTAL_CNT=$(wc -l < $BAM_BED)
echo "total: $TOTAL_CNT"

# rRNA reads
RRNA_CNT=$(bedtools intersect -a $BAM_BED -b $RRNA_BED -u | wc -l)
echo "rRNA: $RRNA_CNT"

# mtRNA reads
MT_CNT=$(bedtools intersect -a $BAM_BED -b $MT_BED -u | wc -l)
echo "mtRNA: $MT_CNT"

# Make exons, introns, genic, intergenic dynamically
EXON_BED=$ANNOTATION/exons.bed
cut -f1-3 $GENOME_BED12 > $EXON_BED
sortBed -i $EXON_BED | mergeBed > $EXON_BED.merged

INTRON_BED=$ANNOTATION/introns.bed
bedtools subtract -a $GENOME_BED12 -b $EXON_BED.merged > $INTRON_BED

GENIC_BED=$ANNOTATION/genic.bed
cat $EXON_BED.merged $INTRON_BED | sortBed | mergeBed > $GENIC_BED

# intergenic
GENOME_SIZE=$ANNOTATION/ChromInfo.txt
INTERGENIC_BED=$ANNOTATION/intergenic.bed
bedtools complement -i $GENIC_BED -g $GENOME_SIZE > $INTERGENIC_BED

# intergenic_notneargene: intergenic at least 5kb from genes
INTERGENIC_NOTNEAR_BED=$ANNOTATION/intergenic_notneargene.bed
bedtools slop -i $GENIC_BED -g $GENOME_SIZE -b 5000 | sortBed | mergeBed \
    | bedtools complement -i - -g $GENOME_SIZE > $INTERGENIC_NOTNEAR_BED

# total_non_rRNA_mt reads
NON_RRNA_MT_BED=$ANNOTATION/non_rRNA_mt.bed
cat $MT_BED $RRNA_BED | sortBed | mergeBed | bedtools complement -i - -g $GENOME_SIZE > $NON_RRNA_MT_BED
NON_RRNA_MT_BAM=$INPUT_BAM.non-rRNA-mt.bam
samtools view -b -L $NON_RRNA_MT_BED $INPUT_BAM > $NON_RRNA_MT_BAM
samtools index $NON_RRNA_MT_BAM
NON_RRNA_MT_CNT=$(samtools view -c $NON_RRNA_MT_BAM)
echo "total_non_rRNA_mt: $NON_RRNA_MT_CNT"

# count reads in all relevant regions
for REGION in intergenic intergenic_notneargene introns exons utr5 utr3 LINE SINE; do
    BED=$ANNOTATION/${REGION}.bed
    CNT=0
    [ -f $BED ] && CNT=$(bedtools intersect -a $BAM_BED -b $BED -u | wc -l)
    echo "$REGION: $CNT"
done
