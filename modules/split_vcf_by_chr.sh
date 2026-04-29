#!/bin/bash
#SBATCH --job-name=split_vcf_by_chr
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

module load BCFtools/1.21

# Paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"

mkdir -p $OUTPUT_DIR

# --- STEP 1: PREP REQUISITES ---
echo "Creating ploidy rules and sex map..."

# A. Create Ploidy Rules (ChrX is ploidy 1 for Males)
RULES="${OUTPUT_DIR}/ploidy_rules.txt"
# Using the full length of ChrX to cover all variants
echo -e "chrX\t1\t156030895\tM\t1" > $RULES

# B. Create Sex Map (Sample <TAB> M/F)
RAW_PASS="${OUTPUT_DIR}/raw_pass_ids.txt"
awk 'NR>1 {print $1}' $PSAM_FILE > $RAW_PASS

SEX_MAP="${OUTPUT_DIR}/sex_map.txt"
# Extract unique subjects and assign M/F
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, (toupper($5)=="MALE" ? "M" : "F")}' | sort -u > ${OUTPUT_DIR}/full_meta.txt
grep -Fwf $RAW_PASS ${OUTPUT_DIR}/full_meta.txt | sed 's/ /\t/' > $SEX_MAP

# Final sample list for subsetting
PASS_SAMPLES="${OUTPUT_DIR}/pass_samples.txt"
cut -f1 $SEX_MAP > $PASS_SAMPLES

# --- STEP 2: CHRX FIX ---
echo "Running fixploidy on chrX..."
# We apply the ploidy fix first, then subset to our 414 samples
bcftools +fixploidy $INPUT_VCF -r chrX -- -s $SEX_MAP -p $RULES | \
bcftools view -S $PASS_SAMPLES --force-samples -Oz -o ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# --- STEP 3: ALL OTHER CHROMOSOMES ---
for chr in {1..22} Y; do
    echo "Processing chr${chr}..."
    bcftools view -r chr${chr} -S $PASS_SAMPLES --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

# --- STEP 4: CLEANUP ---
rm -f ${OUTPUT_DIR}/raw_pass_ids.txt ${OUTPUT_DIR}/full_meta.txt
echo "Success. Pipeline complete."
