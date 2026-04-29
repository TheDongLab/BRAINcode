#!/bin/bash
#SBATCH --job-name=split_vcf_QC_filtered
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

#########################################################################
# This script:
# 1. Extracts the "Pass QC" samples from the PLINK2 filtered dataset.
# 2. Creates a clean male list for those samples.
# 3. Splits the joint VCF by chromosome (1-22, X, Y).
# 4. Forces haploid genotypes for males on ChrX to satisfy TOPMed/eQTL requirements.
#########################################################################

module load BCFtools

# Define paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"

mkdir -p $OUTPUT_DIR

# --- STEP 1: SUBJECT LISTS ---
PASS_SAMPLES="${OUTPUT_DIR}/pass_samples.txt"
awk 'NR>1 {print $1}' $PSAM_FILE > $PASS_SAMPLES

FEMALE_LIST="${OUTPUT_DIR}/females.txt"
MALE_LIST="${OUTPUT_DIR}/males.txt"

# Collapse metadata to unique Subject IDs
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/subject_sex_map.txt
grep -Fwf $PASS_SAMPLES ${OUTPUT_DIR}/subject_sex_map.txt | awk '$2=="female" {print $1}' > $FEMALE_LIST
grep -Fwf $PASS_SAMPLES ${OUTPUT_DIR}/subject_sex_map.txt | awk '$2=="male" {print $1}' > $MALE_LIST

# --- STEP 2: CHRX PROCESS ---
echo "Processing chrX: Isolation Method..."

# 2a. Females: Pure extract (Guaranteed diploid)
bcftools view -r chrX -S $FEMALE_LIST --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/tmp_X_females.vcf.gz
bcftools index -f -t ${OUTPUT_DIR}/tmp_X_females.vcf.gz

# 2b. Males: Force Haploid
# Added the explicit -i filter to satisfy BCFtools 1.21 query requirements
bcftools view -r chrX -S $MALE_LIST --force-samples $INPUT_VCF -Ou | \
bcftools +setGT -- -t q -i 'GT~"." || GT!~"."' -n c:1 | \
bcftools view -Oz -o ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz
bcftools index -f -t ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz

# 2c. Final Merge
bcftools merge --force-samples \
    ${OUTPUT_DIR}/tmp_X_females.vcf.gz \
    ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz \
    -Oz -o ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz
bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# --- STEP 3: ALL OTHER CHRS ---
for chr in {1..22} Y; do
    echo "Processing chr${chr}..."
    bcftools view -r chr${chr} -S $PASS_SAMPLES --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

# --- STEP 4: CLEANUP ---
rm ${OUTPUT_DIR}/tmp_X_females.vcf.gz* ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz*
echo "DONE."
