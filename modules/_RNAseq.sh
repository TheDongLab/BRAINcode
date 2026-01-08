#!/bin/bash
###########################################
# RNAseq pipeline (paired-end)
# Author: Xianjun Dong & Zachery Wolfe
# Date: 1/8/2026
# Version: 2.5
###########################################

set -euo pipefail

echo "[`date`] STEP 0. define paths & inputs"

###########################################
# inputs
###########################################
R1="$1"
R2="$2"
outputdir="$3"

mkdir -p "$outputdir"

CPU=${CPU:-4}
samplename=$(basename "$outputdir")

###########################################
# paths
###########################################
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

GENOME_DIR=/home/zw529/donglab/pipelines/genome
GENOME_FASTA=$GENOME_DIR/hg38.fa
GENOME_STAR_INDEX=$GENOME_DIR
GENOME_GTF=$GENOME_DIR/hg38.ncbiRefSeq.gtf
GENOME_REF_FLAT=$GENOME_DIR/refFlat.txt
GENOME_BED12=$GENOME_DIR/hg38.bed12
CHROM_INFO=$GENOME_DIR/ChromInfo.txt
CHR_M_BED=$GENOME_DIR/chrM.bed
MASK_GTF=""

MAX_HIT=10
strandoption=""

SAMPLE_DIR="$outputdir"
mkdir -p "$SAMPLE_DIR"

###########################################
echo "[`date`] STEP 1. quality check (FastQC)"
###########################################
module load FastQC

[ ! -f "$SAMPLE_DIR/.status.RNAseq.fastqc" ] && \
fastqc -t $CPU --nogroup -o "$SAMPLE_DIR" "$R1" "$R2" && \
touch "$SAMPLE_DIR/.status.RNAseq.fastqc"

###########################################
echo "[`date`] STEP 2. adaptor/quality trimming (Trim Galore)"
###########################################
module load Trim_Galore

[ ! -f "$SAMPLE_DIR/.status.RNAseq.trim" ] && \
trim_galore --paired --cores $CPU --fastqc --output_dir "$SAMPLE_DIR" "$R1" "$R2" && \
touch "$SAMPLE_DIR/.status.RNAseq.trim"

# update FASTQ paths to trimmed files
R1="$SAMPLE_DIR/$(basename "$R1" .fastq.gz)_val_1.fq.gz"
R2="$SAMPLE_DIR/$(basename "$R2" .fastq.gz)_val_2.fq.gz"

###########################################
echo "[`date`] STEP 3. optional kPAL QC"
###########################################
[ ! -f "$SAMPLE_DIR/.status.RNAseq.kpal" ] && \
kpal count "$R1" "$R2" -o "$SAMPLE_DIR/kpal_results" && \
touch "$SAMPLE_DIR/.status.RNAseq.kpal"

###########################################
echo "[`date`] STEP 4. mapping with STAR"
###########################################
module load STAR

STAR_TMP_DIR="$SAMPLE_DIR/STARtmp_${SLURM_JOB_ID:-NA}_${SLURM_ARRAY_TASK_ID:-NA}"

[ ! -f "$SAMPLE_DIR/.status.RNAseq.mapping" ] && \
STAR \
  --runThreadN $CPU \
  --genomeDir "$GENOME_STAR_INDEX" \
  --readFilesIn "$R1" "$R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$SAMPLE_DIR/" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattrRGline ID:$samplename SM:$samplename LB:$samplename PL:ILLUMINA PU:$samplename \
  --outFilterMultimapNmax $MAX_HIT \
  --alignEndsType EndToEnd \
  --chimSegmentMin 10 \
  --chimJunctionOverhangMin 10 \
  --chimOutType Junctions \
  --outSAMstrandField intronMotif \
  --outTmpDir "$STAR_TMP_DIR" && \
touch "$SAMPLE_DIR/.status.RNAseq.mapping"

rm -rf "$STAR_TMP_DIR"

###########################################
echo "[`date`] STEP 4.1. circRNA calling (CIRCexplorer2)"
###########################################
CIRCEXPLORER2=~/donglab/pipelines/modules/rnaseq/bin/CIRCexplorer2
cd "$SAMPLE_DIR"

[ ! -f .status.RNAseq.circRNA ] && \
$CIRCEXPLORER2 parse -t STAR -b Chimeric.out.junction > .CIRCexplorer2_parse.log && \
$CIRCEXPLORER2 annotate \
  -r "$GENOME_REF_FLAT" \
  -g "$GENOME_FASTA" \
  -b back_spliced_junction.bed \
  -o circularRNA_known.txt \
  --low-confidence > .CIRCexplorer2_annotate.log && \
touch .status.RNAseq.circRNA

###########################################
echo "[`date`] STEP 5. post-processing"
###########################################

[ ! -f .status.RNAseq.sam2bam ] && \
samtools index Aligned.sortedByCoord.out.bam && \
touch .status.RNAseq.sam2bam

[ ! -f .status.RNAseq.bam2stat ] && \
echo "$(samtools view -cF 0x100 Aligned.sortedByCoord.out.bam) primary alignments" > Aligned.sortedByCoord.out.bam.stat && \
samtools flagstat Aligned.sortedByCoord.out.bam >> Aligned.sortedByCoord.out.bam.stat && \
touch .status.RNAseq.bam2stat

[ ! -f .status.RNAseq.bam2annotation ] && \
"$SCRIPT_DIR/_bam2annotation.sh" Aligned.sortedByCoord.out.bam > Aligned.sortedByCoord.out.bam.bam2annotation && \
Rscript "$SCRIPT_DIR/_bam2annotation.r" \
  Aligned.sortedByCoord.out.bam.bam2annotation \
  Aligned.sortedByCoord.out.bam.bam2annotation.pdf && \
touch .status.RNAseq.bam2annotation

[ ! -f .status.RNAseq.sam2bw ] && \
total_mapped_reads=$(grep -w total_non_rRNA_mt Aligned.sortedByCoord.out.bam.bam2annotation | cut -f2 -d' ') && \
"$SCRIPT_DIR/bam2bigwig.sh" Aligned.sortedByCoord.out.bam "" "$total_mapped_reads" && \
touch .status.RNAseq.sam2bw

[ ! -f .status.RNAseq.callSNP ] && \
"$SCRIPT_DIR/_callSNP.sh" Aligned.sortedByCoord.out.bam && \
touch .status.RNAseq.callSNP

###########################################
echo "[`date`] STEP 6. transcript assembly / quantification (cufflinks)"
###########################################
module load cufflinks

[ ! -f .status.RNAseq.cufflinks ] && \
cufflinks \
  --no-update-check \
  $strandoption \
  -o "$SAMPLE_DIR" \
  -p $CPU \
  -G "$GENOME_GTF" \
  -M "$MASK_GTF" \
  --compatible-hits-norm \
  --multi-read-correct \
  Aligned.sortedByCoord.out.bam && \
touch .status.RNAseq.cufflinks

###########################################
echo "[`date`] STEP 7. gene counting (htseq-count)"
###########################################

[ ! -f .status.RNAseq.htseqcount ] && \
htseq-count \
  -m intersection-strict \
  -t exon \
  -i gene_id \
  -s no \
  -q \
  -f bam \
  Aligned.sortedByCoord.out.bam \
  "$GENOME_GTF" \
  > hgseqcount.by.gene.tab \
  2> hgseqcount.by.gene.tab.stderr && \
touch .status.RNAseq.htseqcount

echo "[`date`] RNAseq pipeline finished for $samplename"
