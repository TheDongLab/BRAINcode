#!/bin/bash
#SBATCH --job-name=split_vcf_by_chr
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

export LC_ALL=C
module load BCFtools/1.21

# Paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
BASE_OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/chromosome_joint_vcfs"

mkdir -p ${OUTPUT_DIR}

# 1. Generate Sample Lists
# Generate pass_ids directly from the psam file
awk 'NR>1 {print $1}' $PSAM_FILE > ${OUTPUT_DIR}/pass_ids.txt

# Create metadata mappings for Males/Females
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/full_meta.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="male" {print $1}' > ${OUTPUT_DIR}/males.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="female" {print $1}' > ${OUTPUT_DIR}/females.txt

echo "QC samples: $(wc -l < ${OUTPUT_DIR}/pass_ids.txt), Males: $(wc -l < ${OUTPUT_DIR}/males.txt), Females: $(wc -l < ${OUTPUT_DIR}/females.txt)"

# 2. Process Autosomes (1-22) and Y
echo "Processing autosomes and chrY..."
for chr in {1..22} Y; do
    bcftools view -r chr${chr} -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

# 3. Process chrX (Standard Extraction + Diploid Forcing)
echo "Processing chrX with diploid conversion..."

# Extract raw chrX
bcftools view -r chrX -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chrX_raw.vcf.gz

# Apply the Diploid Forcing Logic
bcftools view ${OUTPUT_DIR}/target_ALS_chrX_raw.vcf.gz | \
awk -F'\t' 'BEGIN{OFS="\t"} {
    if(/^#/) {print $0} 
    else {
        for(i=10; i<=NF; i++) {
            if($i == ".") $i="./.";
            else if($i ~ /^\.:/) sub(/^\./, "./.", $i);
            else if($i ~ /^[0-9]+$/) $i=$i"/"$i;
            else if($i ~ /^[0-9]+:/) {
                split($i, a, ":");
                $i = a[1]"/"a[1];
                for(j=2; j<=length(a); j++) $i = $i":"a[j];
            }
        }
        print $0
    }
}' | bcftools view -Oz -o ${OUTPUT_DIR}/target_ALS_chrX_all_diploid.vcf.gz

# Index the final diploid chrX
bcftools index -t ${OUTPUT_DIR}/target_ALS_chrX_all_diploid.vcf.gz

# Clean up raw file to save space
rm ${OUTPUT_DIR}/target_ALS_chrX_raw.vcf.gz

echo "Script complete. Output files are in ${OUTPUT_DIR}"
