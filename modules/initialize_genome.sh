#!/bin/bash
#SBATCH --job-name=init_genome
#SBATCH --output=init_genome.out
#SBATCH --error=init_genome.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --partition=day

module load STAR

# Run STAR in the genome directory
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir . \
     --genomeFastaFiles /home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
     --sjdbGTFfile /home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf \
     --sjdbOverhang 99
