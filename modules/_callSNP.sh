#!/bin/bash
#SBATCH --job-name=callSNP
#SBATCH --output=callSNP.out
#SBATCH --error=callSNP.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G

####################################
## Detect SNP/indels from RNA-seq BAM/SAM
## Updated to use genome folder: /home/zw529/donglab/pipelines/genome
####################################

if [ $# -ne 1 ]; then
    echo "Usage: `basename $0` <input.sam or .bam>"
    exit
fi

input_sam=$1
pipeline_path=/home/zw529/donglab/pipelines/modules/rnaseq
export PATH=$pipeline_path:$PATH

GENOME_DIR=/home/zw529/donglab/pipelines/genome

dbSNP=$GENOME_DIR/dbSNP.bed
rRNA=$GENOME_DIR/rRNA.bed
LINE=$GENOME_DIR/LINE.bed
exons=$GENOME_DIR/exons.bed
introns=$GENOME_DIR/introns.bed

# convert BAM to SAM if necessary
if [[ $input_sam == *.bam ]]; then
    samtools view -h -o ${input_sam%.bam}.sam $input_sam
    input_sam=${input_sam%.bam}.sam
fi

# Split large files
size=$(ls -sH $input_sam | cut -f1 -d' ')
echo "File size is $size KB"

if [ "$size" -gt 10000000 ]; then
    split -b 1000M $input_sam tmp_sampiece_
    > .paraFile
    for i in tmp_sampiece_*; do
        echo "awk -f $pipeline_path/modules/_sam2variation.awk $i > tmp_snp.$i" >> .paraFile
    done
    rm -f .paraFile.completed
    ParaFly -c .paraFile -CPU 8

    cat tmp_snp* | sort -u > ${input_sam/sam/snp}
    rm tmp_snp* tmp_sampiece_*
else
    awk -f $pipeline_path/modules/_sam2variation.awk $input_sam | sort -u > ${input_sam/sam/snp}
fi

# SNP position histogram
textHistogram -col=4 -maxBinCount=100 ${input_sam/sam/snp} > ${input_sam/sam/snp}.relpos.hist

# get read length
readslength=$(grep -v -m 10 "^@" $input_sam | awk '{print length($10)}' | sort -nr | head -n1)

# filter SNPs near read ends
awk -v RL=$readslength '$4>=10 && $4<=(RL-10)' ${input_sam/sam/snp} \
    | cut -f1-2 \
    | sed 's/:/\t/;s/-/\t/' \
    | bedtools groupby -g 1-4 -c 4 -o count \
    | sort -k5,5nr > ${input_sam/sam/snp.depth}

textHistogram -col=5 -minVal=10 -maxBinCount=100 -binSize=10 ${input_sam/sam/snp.depth} > ${input_sam/sam/snp.depth.hist

# retain SNPs with depth >=16
awk '$5>15' ${input_sam/sam/snp}.depth | sortBed > ${input_sam/sam/snp}.depth_gt_15

# annotate SNPs
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $dbSNP   -sorted -wao | cut -f1-5,9 | groupBy -g 1-5 -c 6 -o collapse > ${input_sam/sam/snp}_dbSNP
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $rRNA    -sorted -wao | cut -f1-5,9 | groupBy -g 1-5 -c 6 -o collapse > ${input_sam/sam/snp}_rRNA
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $LINE    -sorted -wao | cut -f1-5,9 | groupBy -g 1-5 -c 6 -o collapse > ${input_sam/sam/snp}_LINE
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $exons   -sorted -wao | cut -f1-5,9 | uniq | groupBy -g 1-5 -c 6 -o collapse > ${input_sam/sam/snp}_exon
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $introns -sorted -wao | cut -f1-5,9 | uniq | groupBy -g 1-5 -c 6 -o collapse > ${input_sam/sam/snp}_intron

# combine annotations in order: dbSNP, exon, intron, LINE, rRNA
paste ${input_sam/sam/snp}_* | awk '{OFS="\t"; printf "%s\t%s\t%s\t%s\t%s",$1,$2,$3,$4,$5; for(i=6;i<=NF;i+=6) printf "\t%s",$i; printf "\n"}' > ${input_sam/sam/snp}.annotation

# cleanup
rm ${input_sam/sam/snp}_*

echo "SNP annotation complete: ${input_sam/sam/snp}.annotation"
