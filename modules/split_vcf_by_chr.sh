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

# --- STEP 1: GENERATE CLEAN SUBJECT LISTS ---
# We extract the 414 subjects that passed PLINK QC
PASS_SAMPLES="${OUTPUT_DIR}/pass_samples.txt"
awk 'NR>1 {print $1}' $PSAM_FILE > $PASS_SAMPLES

# Fix the metadata pull: 
# 1. tr -d '\r' removes Windows line endings
# 2. awk uses Column 2 (SubjectID) and Column 5 (Sex)
# 3. sort -u collapses the thousands of tissue samples into unique subject entries
echo "Collapsing metadata to subject-level..."
FEMALE_LIST="${OUTPUT_DIR}/females.txt"
MALE_LIST="${OUTPUT_DIR}/males.txt"

tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/subject_sex_map.txt

# Now intersect that unique map with your 414 QC-pass list
grep -Fwf $PASS_SAMPLES ${OUTPUT_DIR}/subject_sex_map.txt | awk '$2=="female" {print $1}' > $FEMALE_LIST
grep -Fwf $PASS_SAMPLES ${OUTPUT_DIR}/subject_sex_map.txt | awk '$2=="male" {print $1}' > $MALE_LIST

echo "Stats: Pass QC: $(wc -l < $PASS_SAMPLES) | Unique Males: $(wc -l < $MALE_LIST) | Unique Females: $(wc -l < $FEMALE_LIST)"

# --- STEP 2: PROCESS CHROMOSOME X (The Sex-Ploidy Fix) ---
echo "Processing chrX..."

# Extract Females
bcftools view -r chrX -S $FEMALE_LIST --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/tmp_X_females.vcf.gz
bcftools index -t ${OUTPUT_DIR}/tmp_X_females.vcf.gz

# Extract Males and Haploidize
bcftools view -r chrX -S $MALE_LIST --force-samples $INPUT_VCF -Ou | \
bcftools +setGT -- -t q -n h | \
bcftools view -Oz -o ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz
bcftools index -t ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz

# Merge
bcftools merge --force-samples \
    ${OUTPUT_DIR}/tmp_X_females.vcf.gz \
    ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz \
    -Oz -o ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# --- STEP 3: OTHERS ---
for chr in {1..22} Y; do
    echo "Processing chr${chr}..."
    bcftools view -r chr${chr} -S $PASS_SAMPLES --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

# Cleanup
rm ${OUTPUT_DIR}/tmp_X_females.vcf.gz* ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz*
