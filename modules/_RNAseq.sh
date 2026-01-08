#!/bin/bash
###########################################
# RNAseq pipeline (paired-end)
# Author: Xianjun Dong & Zachery Wolfe
# Date: 1/8/2026
# Version: 2.2
###########################################

echo "[`date`] STEP 0. define paths & inputs"

# Third argument: output folder
outputdir="$3"
mkdir -p "$outputdir"

# Number of threads
CPU=${CPU:-4}

# Sample name inferred from folder name
samplename=$(basename "$outputdir")

# Input FASTQ (must be full paths)
R1="$1"
R2="$2"

# Paths
SCRIPT_DIR=$(dirname "$0")
GENOME_DIR=/home/zw529/donglab/pipelines/genome
GENOME_FASTA=$GENOME_DIR/hg38.fa
GENOME_STAR_INDEX=$GENOME_DIR
GENOME_GTF=$GENOME_DIR/hg38.ncbiRefSeq.gtf
GENOME_REF_FLAT=$GENOME_DIR/refFlat.txt
GENOME_BED12=$GENOME_DIR/hg38.bed12
CHROM_INFO=$GENOME_DIR/ChromInfo.txt
CHR_M_BED=$GENOME_DIR/chrM.bed
MASK_GTF=""  # optional

MAX_HIT=10
strandoption=""

mkdir -p "$outputdir/$samplename"

###########################################
echo "[`date`] STEP 1. quality check (FastQC)"
###########################################
module load FastQC

[ ! -f "$outputdir/$samplename/.status.RNAseq.fastqc" ] && \
fastqc -t $CPU --nogroup -o "$outputdir/$samplename" "$R1" "$R2" && \
touch "$outputdir/$samplename/.status.RNAseq.fastqc"

###########################################
echo "[`date`] STEP 2. adaptor/quality trimming (Trim Galore)"
###########################################
module load Trim_Galore

[ ! -f "$outputdir/$samplename/.status.RNAseq.trim" ] && \
trim_galore --paired --cores $CPU --fastqc --output_dir "$outputdir/$samplename" "$R1" "$R2" && \
touch "$outputdir/$samplename/.status.RNAseq.trim"

# Update R1/R2 to trimmed outputs
R1="$outputdir/$samplename/$(basename $R1 .fastq.gz)_val_1.fq.gz"
R2="$outputdir/$samplename/$(basename $R2 .fastq.gz)_val_2.fq.gz"

###########################################
echo "[`date`] STEP 3. optional kPAL QC / library complexity"
###########################################
[ ! -f "$outputdir/$samplename/.status.RNAseq.kpal" ] && \
kpal count "$R1" "$R2" -o "$outputdir/$samplename/kpal_results" && \
touch "$outputdir/$samplename/.status.RNAseq.kpal"

###########################################
echo "[`date`] STEP 4. mapping with STAR"
###########################################
module load STAR

[ ! -f "$outputdir/$samplename/.status.RNAseq.mapping" ] && \
STAR --runThreadN $CPU \
     --genomeDir $GENOME_STAR_INDEX \
     --readFilesIn "$R1" "$R2" \
     --readFilesCommand zcat \
     --outFileNamePrefix "$outputdir/$samplename/" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattrRGline ID:$samplename SM:$samplename LB:$samplename PL:ILLUMINA PU:$samplename \
     --outFilterMultimapNmax $MAX_HIT \
     --alignEndsType EndToEnd \
     --chimSegmentMin 10 \
     --chimJunctionOverhangMin 10 \
     --chimOutType Junctions \
     --outSAMstrandField intronMotif && \
touch "$outputdir/$samplename/.status.RNAseq.mapping"

###########################################
echo "[`date`] STEP 4.1. circRNA calling with CIRCexplorer2"
###########################################
module load CIRCexplorer2

cd "$outputdir/$samplename"

[ ! -f .status.RNAseq.circRNA ] && \
CIRCexplorer2 parse -t STAR -b Chimeric.out.junction > .CIRCexplorer2_parse.log && \
CIRCexplorer2 annotate -r $GENOME_REF_FLAT -g $GENOME_FASTA -b back_spliced_junction.bed \
     -o circularRNA_known.txt --low-confidence > .CIRCexplorer2_annotate.log && \
touch .status.RNAseq.circRNA

[ ! -f .status.RNAseq.circRNA.circexplore.denovo ] && \
CIRCexplorer2 assemble -r $GENOME_REF_FLAT -m . -o assemble -p $CPU --remove-rRNA > .CIRCexplorer2_assemble.log && \
CIRCexplorer2 denovo -r $GENOME_REF_FLAT --as=AS --abs=ABS -g $GENOME_FASTA \
     -b back_spliced_junction.bed -d assemble -m . -o denovo_circ > .CIRCexplorer2_denovo.log && \
touch .status.RNAseq.circRNA.circexplore.denovo

###########################################
echo "[`date`] STEP 5. post-processing, indexing, stats, bigwig"
###########################################

cd "$outputdir/$samplename"

# index STAR BAM
[ ! -f .status.RNAseq.sam2bam ] && \
samtools index Aligned.sortedByCoord.out.bam && \
touch .status.RNAseq.sam2bam

# alignment stats
[ ! -f .status.RNAseq.bam2stat ] && \
echo "$(samtools view -cF 0x100 Aligned.sortedByCoord.out.bam) primary alignments" > Aligned.sortedByCoord.out.bam.stat && \
samtools flagstat Aligned.sortedByCoord.out.bam >> Aligned.sortedByCoord.out.bam.stat && \
touch .status.RNAseq.bam2stat

# bam2annotation
[ ! -f .status.RNAseq.bam2annotation ] && \
$SCRIPT_DIR/_bam2annotation.sh Aligned.sortedByCoord.out.bam > Aligned.sortedByCoord.out.bam.bam2annotation && \
Rscript $SCRIPT_DIR/_bam2annotation.r Aligned.sortedByCoord.out.bam.bam2annotation \
     Aligned.sortedByCoord.out.bam.bam2annotation.pdf && \
touch .status.RNAseq.bam2annotation

# bigwig
[ ! -f .status.RNAseq.sam2bw ] && \
total_mapped_reads=$(grep -w total_non_rRNA_mt Aligned.sortedByCoord.out.bam.bam2annotation | cut -f2 -d' ') && \
$SCRIPT_DIR/bam2bigwig.sh Aligned.sortedByCoord.out.bam "" $total_mapped_reads && \
touch .status.RNAseq.sam2bw

# SNP calling
[ ! -f .status.RNAseq.callSNP ] && \
$SCRIPT_DIR/_callSNP.sh Aligned.sortedByCoord.out.bam && \
touch .status.RNAseq.callSNP

# assembly/quantification with cufflinks
[ ! -f .status.RNAseq.cufflinks ] && \
cufflinks --no-update-check $strandoption -o ./ -p $CPU -G $GENOME_GTF -M $MASK_GTF \
     --compatible-hits-norm --multi-read-correct Aligned.sortedByCoord.out.bam && \
touch .status.RNAseq.cufflinks

# gene counting
[ ! -f .status.RNAseq.htseqcount ] && \
htseq-count -m intersection-strict -t exon -i gene_id -s no -q -f bam \
     Aligned.sortedByCoord.out.bam $GENOME_GTF > hgseqcount.by.gene.tab 2> hgseqcount.by.gene.tab.stderr && \
touch .status.RNAseq.htseqcount

echo "[`date`] RNAseq pipeline finished for $samplename"
