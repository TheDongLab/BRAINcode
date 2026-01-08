#!/bin/bash
#SBATCH --job-name=callSNP
#SBATCH --output=callSNP.out
#SBATCH --error=callSNP.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G

###########################################
# Detect SNP/indels from RNA-seq BAM/SAM
# Updated for local genome and rnaseq modules
###########################################

set -e

if [ $# -ne 1 ]; then
    echo "Usage: $(basename $0) <input.sam or .bam>"
    exit 1
fi

input_sam="$1"
pipeline_path=/home/zw529/donglab/pipelines/scripts/rnaseq
module_path=/home/zw529/donglab/pipelines/scripts/rnaseq

GENOME_DIR=/home/zw529/donglab/pipelines/genome
dbSNP=$GENOME_DIR/dbSNP.bed
rRNA=$GENOME_DIR/rRNA.bed
LINE=$GENOME_DIR/LINE.bed
exons=$GENOME_DIR/exons.bed
introns=$GENOME_DIR/introns.bed

# convert BAM to SAM if necessary
if [[ "$input_sam" == *.bam ]]; then
    samtools view -h -o "${input_sam%.bam}.sam" "$input_sam"
    input_sam="${input_sam%.bam}.sam"
fi

# Split large files if >10 GB
size_kb=$(ls -sH "$input_sam" | awk '{print $1}')
echo "File size is $size_kb KB"

cd "$(dirname "$input_sam")"
bam_base=$(basename "$input_sam" .sam)

if [ "$size_kb" -gt 10000000 ]; then
    split -b 1000M "$input_sam" tmp_sampiece_
    > .paraFile
    for i in tmp_sampiece_*; do
        echo "awk -f $module_path/_sam2variation.awk $i > tmp_snp.$i" >> .paraFile
    done
    rm -f .paraFile.completed
    ParaFly -c .paraFile -CPU 8
    cat tmp_snp* | sort -u > "${bam_base}.snp"
    rm -f tmp_snp* tmp_sampiece_*
else
    awk -f "$module_path/_sam2variation.awk" "$input_sam" | sort -u > "${bam_base}.snp"
fi

# SNP position histogram
textHistogram -col=4 -maxBinCount=100 "${bam_base}.snp" > "${bam_base}.snp.relpos.hist"

# get read length
readslength=$(awk 'BEGIN{max=0} !/^@/{len=length($10); if(len>max) max=len} END{print max}' "$input_sam")

# filter SNPs near read ends
awk -v RL=$readslength '$4>=10 && $4<=(RL-10)' "${bam_base}.snp" \
    | cut -f1-2 \
    | sed 's/:/\t/;s/-/\t/' \
    | bedtools groupby -g 1-4 -c 4 -o count \
    | sort -k5,5nr > "${bam_base}.snp.depth"

textHistogram -col=5 -minVal=10 -maxBinCount=100 -binSize=10 "${bam_base}.snp.depth" > "${bam_base}.snp.depth.hist"

# retain SNPs with depth >=16
awk '$5>15' "${bam_base}.snp.depth" | sortBed > "${bam_base}.snp.depth_gt_15"

# annotate SNPs
intersectBed -a "${bam_base}.snp.depth_gt_15" -b "$dbSNP"   -sorted -wao | cut -f1-5,9 | groupBy -g 1-5 -c 6 -o collapse > "${bam_base}.snp_dbSNP"
intersectBed -a "${bam_base}.snp.depth_gt_15" -b "$rRNA"    -sorted -wao | cut -f1-5,9 | groupBy -g 1-5 -c 6 -o collapse > "${bam_base}.snp_rRNA"
intersectBed -a "${bam_base}.snp.depth_gt_15" -b "$LINE"    -sorted -wao | cut -f1-5,9 | groupBy -g 1-5 -c 6 -o collapse > "${bam_base}.snp_LINE"
intersectBed -a "${bam_base}.snp.depth_gt_15" -b "$exons"   -sorted -wao | cut -f1-5,9 | uniq | groupBy -g 1-5 -c 6 -o collapse > "${bam_base}.snp_exon"
intersectBed -a "${bam_base}.snp.depth_gt_15" -b "$introns" -sorted -wao | cut -f1-5,9 | uniq | groupBy -g 1-5 -c 6 -o collapse > "${bam_base}.snp_intron"

# combine annotations: dbSNP, exon, intron, LINE, rRNA
paste "${bam_base}.snp_"* | awk '{
    OFS="\t";
    printf "%s\t%s\t%s\t%s\t%s",$1,$2,$3,$4,$5;
    for(i=6;i<=NF;i+=6) printf "\t%s",$i;
    printf "\n"
}' > "${bam_base}.snp.annotation"

# cleanup intermediate annotation files
rm -f "${bam_base}.snp_"*

echo "SNP annotation complete: ${bam_base}.snp.annotation"
