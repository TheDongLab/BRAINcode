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

# --- STEP 1: CREATE SAMPLE LISTS ---
PASS_SAMPLES="${OUTPUT_DIR}/pass_samples.txt"
awk 'NR>1 {print $1}' $PSAM_FILE > $PASS_SAMPLES

FEMALE_LIST="${OUTPUT_DIR}/females.txt"
MALE_LIST="${OUTPUT_DIR}/males.txt"

# Clean metadata and create sex map
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/subject_sex_map.txt

grep -Fwf $PASS_SAMPLES ${OUTPUT_DIR}/subject_sex_map.txt | awk '$2=="female" {print $1}' > $FEMALE_LIST
grep -Fwf $PASS_SAMPLES ${OUTPUT_DIR}/subject_sex_map.txt | awk '$2=="male" {print $1}' > $MALE_LIST

echo "Females: $(wc -l < $FEMALE_LIST)"
echo "Males: $(wc -l < $MALE_LIST)"

# --- STEP 2: PROCESS CHRX ---
echo "Processing chrX..."

# 2a. Extract females (keep as-is, diploid)
bcftools view -r chrX -S $FEMALE_LIST --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/tmp_X_females.vcf.gz
bcftools index -t ${OUTPUT_DIR}/tmp_X_females.vcf.gz

# 2b. Extract males as diploid (uncompressed intermediate)
bcftools view -r chrX -S $MALE_LIST --force-samples $INPUT_VCF -Ou -o ${OUTPUT_DIR}/tmp_X_males_diploid.vcf

# 2c. Convert males to haploid using bcftools +setGT
# -t q -i 'GT!="."': apply to all non-missing genotypes
# -n c:1: convert to ploidy 1 (haploid)
bcftools +setGT -Oz -o ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz ${OUTPUT_DIR}/tmp_X_males_diploid.vcf -- -t q -i 'GT!="."' -n c:1

bcftools index -t ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz

# Verify the conversion worked
echo "Female ploidy check (should be diploid like 0/0, 0/1, 1/1):"
bcftools query -f '[%GT\n]' ${OUTPUT_DIR}/tmp_X_females.vcf.gz | head -5

echo "Male ploidy check (should be haploid like 0, 1):"
bcftools query -f '[%GT\n]' ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz | head -5

# 2d. Merge females and haploid males into single VCF
bcftools merge --force-samples \
    ${OUTPUT_DIR}/tmp_X_females.vcf.gz \
    ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz \
    -Oz -o ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# --- STEP 3: ALL OTHER CHROMOSOMES ---
for chr in {1..22} Y; do
    echo "Processing chr${chr}..."
    bcftools view -r chr${chr} -S $PASS_SAMPLES --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

# --- STEP 4: CLEANUP ---
rm -f ${OUTPUT_DIR}/tmp_X_females.vcf.gz*
rm -f ${OUTPUT_DIR}/tmp_X_males_diploid.vcf*
rm -f ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz*

echo "Successfully created mixed-ploidy chrX VCF!"
echo "chrX output: ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz"
