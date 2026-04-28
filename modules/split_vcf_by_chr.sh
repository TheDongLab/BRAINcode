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

# --- STEP 1: GENERATE SAMPLE LISTS ---
echo "Generating sample lists..."
PASS_SAMPLES="${OUTPUT_DIR}/pass_samples.txt"
# Extract the 413 samples from your filtered PLINK results
awk 'NR>1 {print $1}' $PSAM_FILE > $PASS_SAMPLES

# Define Females and Males specifically from the PASS list
FEMALE_LIST="${OUTPUT_DIR}/females.txt"
MALE_LIST="${OUTPUT_DIR}/males.txt"

tr -d '\r' < $METADATA | awk -F',' 'NR>1 && (tolower($5) == "female") {print $2}' | grep -Fwf $PASS_SAMPLES > $FEMALE_LIST
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && (tolower($5) == "male") {print $2}' | grep -Fwf $PASS_SAMPLES > $MALE_LIST

echo "Stats: Total Passing: $(wc -l < $PASS_SAMPLES) | Males: $(wc -l < $MALE_LIST) | Females: $(wc -l < $FEMALE_LIST)"

# --- STEP 2: PROCESS CHROMOSOME X (The Sex-Ploidy Fix) ---
echo "Processing chrX: Splitting by sex for ploidy fix..."

# 2a. Extract Females (Leave as diploid)
bcftools view -r chrX -S $FEMALE_LIST --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/tmp_X_females.vcf.gz
bcftools index -t ${OUTPUT_DIR}/tmp_X_females.vcf.gz

# 2b. Extract Males and Force Haploid
# +setGT -n h (sets genotypes to haploid)
bcftools view -r chrX -S $MALE_LIST --force-samples $INPUT_VCF -Ou | \
bcftools +setGT -- -t q -n h | \
bcftools view -Oz -o ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz
bcftools index -t ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz

# 2c. Merge them back together
echo "Merging chrX back into a single file..."
bcftools merge --force-samples \
    ${OUTPUT_DIR}/tmp_X_females.vcf.gz \
    ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz \
    -Oz -o ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# --- STEP 3: PROCESS ALL OTHER CHROMOSOMES ---
for chr in {1..22} Y; do
    echo "Processing chr${chr}..."
    bcftools view -r chr${chr} -S $PASS_SAMPLES --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

# --- STEP 4: CLEANUP ---
echo "Cleaning up temporary files..."
rm ${OUTPUT_DIR}/tmp_X_females.vcf.gz* ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz*

echo "PIPELINE COMPLETE."
echo "Final files are in: $OUTPUT_DIR"
