###########################################
# bash script for running paired-end RNAseq
# author: Xianjun Dong and Zachery Wolfe (modified)
# email: xdong@rics.bwh.harvard.edu
# date: 1/8/2026
# version: 2.1
# Note: call this script in the folder containing raw FASTQ files (R1.fastq.gz / R2.fastq.gz)
###########################################
#!/bin/bash

###########################################
echo "["`date`"] STEP 0. define paths & inputs"
###########################################

# Paths to binaries and pipeline scripts
BIN_DIR=~/donglab/pipelines/modules/rnaseq/bin
SCRIPT_DIR=~/donglab/pipelines/scripts/rnaseq

outputdir="$3"      # third argument from run_RNAseq.txt
samplename=$(basename "$R1" | sed 's/_1\.fastq\.gz//')  # extract sample name from FASTQ1
CPU=8               # or pass as arg

# Genome references (using initialized STAR hg38)
GENOME_DIR=/home/zw529/donglab/pipelines/genome
GENOME_FASTA=$GENOME_DIR/hg38.fa
GENOME_STAR_INDEX=$GENOME_DIR
GENOME_GTF=$GENOME_DIR/hg38.ncbiRefSeq.gtf
GENOME_REF_FLAT=$GENOME_DIR/refFlat.txt
GENOME_BED12=""

# Optional mask GTF for cufflinks
MASK_GTF=""

# STAR mapping options
MAX_HIT=10
strandoption=""  # set "--library-type fr-firststrand" for stranded data

# Input FASTQ (assumes running inside folder containing R1/R2)
R1=$(ls *_1.fastq.gz | head -n1)
R2=$(ls *_2.fastq.gz | head -n1)

# Sample and output names
samplename=${PWD##*/}        # folder name as sample name
inputdir=$PWD
outputdir=$inputdir/../processed
modulename="RNAseq"

# Make output folder if it doesn't exist
mkdir -p $outputdir/$samplename

###########################################
echo "["`date`"] STEP 1. quality check (FastQC)"
###########################################

module load FastQC
[ ! -f $outputdir/$samplename/.status.$modulename.fastqc ] && \
fastqc -t $CPU -o $outputdir/$samplename $R1 $R2 && \
touch $outputdir/$samplename/.status.$modulename.fastqc

###########################################
echo "["`date`"] STEP 2. adaptor/quality trimming (Trim Galore)"
###########################################

module load TrimGalore

[ ! -f $outputdir/$samplename/.status.$modulename.trim ] && \
trim_galore --paired --cores $CPU --fastqc --output_dir $outputdir/$samplename $R1 $R2 && \
touch $outputdir/$samplename/.status.$modulename.trim

# Update R1/R2 to trimmed files
R1=$outputdir/$samplename/$(basename $R1 .fastq.gz)_val_1.fq.gz
R2=$outputdir/$samplename/$(basename $R2 .fastq.gz)_val_2.fq.gz

###########################################
echo "["`date`"] STEP 3. optional kPAL QC / library complexity"
###########################################

[ ! -f $outputdir/$samplename/.status.$modulename.kpal ] && \
kpal -f $R1 $R2 -o $outputdir/$samplename/kpal_results && \
touch $outputdir/$samplename/.status.$modulename.kpal

###########################################
echo "["`date`"] STEP 4. mapping with STAR"
###########################################

module load STAR

[ ! -f $outputdir/$samplename/.status.$modulename.mapping ] && \
STAR --runThreadN $CPU \
     --genomeDir $GENOME_STAR_INDEX \
     --readFilesIn $R1 $R2 \
     --readFilesCommand zcat \
     --outFileNamePrefix $outputdir/$samplename/ \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattrRGline ID:$samplename SM:$samplename LB:$samplename PL:ILLUMINA PU:$samplename \
     --outFilterMultimapNmax $MAX_HIT \
     --alignEndsType EndToEnd \
     --chimSegmentMin 10 \
     --chimJunctionOverhangMin 10 \
     --chimOutType Junctions \
     --outSAMstrandField intronMotif && \
touch $outputdir/$samplename/.status.$modulename.mapping

###########################################
echo "["`date`"] STEP 4.1. circRNA calling with CIRCexplorer2"
###########################################

module load CIRCexplorer2

cd $outputdir/$samplename

[ ! -f .status.$modulename.circRNA ] && \
CIRCexplorer2 parse -t STAR \
     -b Chimeric.out.junction \
     > .CIRCexplorer2_parse.log && \
CIRCexplorer2 annotate -r $GENOME_REF_FLAT \
     -g $GENOME_FASTA \
     -b back_spliced_junction.bed \
     -o circularRNA_known.txt --low-confidence \
     > .CIRCexplorer2_annotate.log && \
touch .status.$modulename.circRNA

[ ! -f .status.$modulename.circRNA.circexplore.denovo ] && \
CIRCexplorer2 assemble -r $GENOME_REF_FLAT -m . -o assemble -p $CPU --remove-rRNA > .CIRCexplorer_assemble.log && \
CIRCexplorer2 denovo -r $GENOME_REF_FLAT --as=AS --abs=ABS \
     -g $GENOME_FASTA \
     -b back_spliced_junction.bed -d assemble -m . -o denovo_circ \
     > .CIRCexplorer_denovo.log && \
touch .status.$modulename.circRNA.circexplore.denovo

###########################################
echo "["`date`"] STEP 5. post-processing, indexing, stats, bigwig"
###########################################

cd $outputdir/$samplename

# index STAR BAM
[ ! -f .status.$modulename.sam2bam ] && \
samtools index Aligned.sortedByCoord.out.bam && \
touch .status.$modulename.sam2bam

# alignment stats
[ ! -f .status.$modulename.bam2stat ] && \
echo `samtools view -cF 0x100 Aligned.sortedByCoord.out.bam` "primary alignments" > Aligned.sortedByCoord.out.bam.stat && \
samtools flagstat Aligned.sortedByCoord.out.bam >> Aligned.sortedByCoord.out.bam.stat && \
touch .status.$modulename.bam2stat

# bam2annotation
[ ! -f .status.$modulename.bam2annotation ] && \
([ -f Aligned.sortedByCoord.out.bam.bam2annotation ] || $SCRIPT_DIR/_bam2annotation.sh Aligned.sortedByCoord.out.bam > Aligned.sortedByCoord.out.bam.bam2annotation) && \
Rscript $SCRIPT_DIR/_bam2annotation.r Aligned.sortedByCoord.out.bam.bam2annotation Aligned.sortedByCoord.out.bam.bam2annotation.pdf && \
touch .status.$modulename.bam2annotation

# bigwig for UCSC
[ ! -f .status.$modulename.sam2bw ] && \
total_mapped_reads=`grep -w total_non_rRNA_mt Aligned.sortedByCoord.out.bam.bam2annotation | cut -f2 -d' '` && \
bam2bigwig.sh Aligned.sortedByCoord.out.bam "" $total_mapped_reads && \
touch .status.$modulename.sam2bw

# SNP calling
[ ! -f .status.$modulename.callSNP ] && \
$SCRIPT_DIR/_callSNP.sh Aligned.sortedByCoord.out.bam && \
touch .status.$modulename.callSNP

# assembly/quantification with cufflinks and htseq-count
[ ! -f .status.$modulename.cufflinks ] && \
cufflinks --no-update-check $strandoption -o ./ -p $CPU -G $GENOME_GTF -M $MASK_GTF --compatible-hits-norm --multi-read-correct Aligned.sortedByCoord.out.bam && \
touch .status.$modulename.cufflinks

[ ! -f .status.$modulename.htseqcount ] && \
htseq-count -m intersection-strict -t exon -i gene_id -s no -q -f bam Aligned.sortedByCoord.out.bam $GENOME_GTF > hgseqcount.by.gene.tab 2> hgseqcount.by.gene.tab.stderr && \
touch .status.$modulename.htseqcount
