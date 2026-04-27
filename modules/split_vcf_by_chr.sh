#!/bin/bash
#SBATCH --job-name=split_vcf
#SBATCH --output=split_vcf.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00

###########################################
# Batch script to split joint .vcf file by chromosome. This step is necessary to perform imputation (post-QC) on the TOPMed Server
###########################################

module load BCFtools

# Define paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/Cerebellum/eQTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through chromosomes 1-22, X, and Y
for chr in {1..22} X Y; do
    echo "Processing chromosome: chr${chr}..."
    
    # Extract specific chromosome, compress, and save to the new folder
    # Note: Use "${chr}" or "chr${chr}" depending on your VCF's naming convention
    bcftools view -r chr${chr} $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    
    # Index the output file (required for the imputation server)
    bcftools index -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

echo "Splitting complete. Files are located in: $OUTPUT_DIR"
