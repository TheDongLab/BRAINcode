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

# Define paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"

mkdir -p $OUTPUT_DIR

# --- STEP 1: SAMPLE LISTS & PLOIDY MAP ---
echo "Generating sample lists and ploidy map..."

# Extract IDs from PLINK pass list
RAW_PASS="${OUTPUT_DIR}/raw_pass_ids.txt"
awk 'NR>1 {print $1}' $PSAM_FILE > $RAW_PASS

# Create BCFtools sex map: SampleID [TAB] Sex (1 for male, 2 for female)
# We intersect the metadata with the PLINK pass list to ensure perfect sample consistency
PLOIDY_MAP="${OUTPUT_DIR}/ploidy_map.txt"
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, (tolower($5)=="male" ? "1" : "2")}' | sort -u > ${OUTPUT_DIR}/full_sex_map.txt
grep -Fwf $RAW_PASS ${OUTPUT_DIR}/full_sex_map.txt | sed 's/ /\t/' > $PLOIDY_MAP

# Update PASS_SAMPLES to only include those in the ploidy map (fixes the 411/414 mismatch)
PASS_SAMPLES="${OUTPUT_DIR}/pass_samples.txt"
cut -f1 $PLOIDY_MAP > $PASS_SAMPLES

echo "Final Subject Count: $(wc -l < $PASS_SAMPLES)"

# --- STEP 2: PROCESS CHRX (THE FIX) ---
echo "Processing chrX with fixploidy..."

# We subset to the 414/411 samples and immediately apply fixploidy
# This forces Sex 1 (Male) to haploid and Sex 2 (Female) to diploid
bcftools view -r chrX -S $PASS_SAMPLES --force-samples $INPUT_VCF -Ou | \
bcftools +fixploidy -- -p $PLOIDY_MAP | \
bcftools view -Oz -o ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# --- STEP 3: ALL OTHER CHROMOSOMES ---
for chr in {1..22} Y; do
    echo "Processing chr${chr}..."
    bcftools view -r chr${chr} -S $PASS_SAMPLES --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

# --- STEP 4: CLEANUP ---
rm -f ${OUTPUT_DIR}/raw_pass_ids.txt ${OUTPUT_DIR}/full_sex_map.txt
echo "Success."
