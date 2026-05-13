#!/bin/bash
#SBATCH --job-name=TargetALS_ChrX_Full_QC
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/chrX_comprehensive_%j.out

# --- PATHS ---
VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
PSAM_IN="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"
REF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"

# Your source lists
USER_DATA_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
MALES_SRC="${USER_DATA_DIR}/males.txt"
FEMALES_SRC="${USER_DATA_DIR}/females.txt"

module load BCFtools/1.21
module load PLINK/1.90-beta6.10

# --- STEP 1: PREP SEX UPDATE & AUDIT ---

# 1. Create the cohort keep list from .psam
awk 'NR>1 {print $1, $2}' $PSAM_IN > ${OUT_DIR}/cohort_keep_plink.txt

# 2. Create the 3-column sex update file (FID IID SEX) for PLINK
# Using --double-id logic: duplicating the sample ID for both FID and IID
awk '{print $1, $1, "1"}' $MALES_SRC > ${OUT_DIR}/sex_update.txt
awk '{print $1, $1, "2"}' $FEMALES_SRC >> ${OUT_DIR}/sex_update.txt

# 3. Run PLINK sex check with split-x for hg38
# This validates your lists against the actual X-chromosome heterozygosity
plink --vcf $VCF_IN \
      --double-id \
      --keep ${OUT_DIR}/cohort_keep_plink.txt \
      --update-sex ${OUT_DIR}/sex_update.txt \
      --split-x hg38 \
      --check-sex \
      --allow-extra-chr \
      --out ${OUT_DIR}/genetic_sex_audit

# 4. Extract IDs for subsetting, skipping the "IID" header row
# We trust the genetic F-statistic here to ensure pure cohorts
awk 'NR>1 && $6 > 0.8 {print $2}' ${OUT_DIR}/genetic_sex_audit.sexcheck > ${OUT_DIR}/final_males.txt
awk 'NR>1 && $6 < 0.2 {print $2}' ${OUT_DIR}/genetic_sex_audit.sexcheck > ${OUT_DIR}/final_females.txt

# --- STEP 2: NORMALIZATION & STRATIFIED SPLIT ---

for SEX in males females; do
    echo "Processing ${SEX^^}..."
    
    ID_FILE="${OUT_DIR}/final_${SEX}.txt"
    
    if [ ! -s "$ID_FILE" ]; then
        echo "Warning: No samples passed audit for ${SEX}. Check .sexcheck file."
        continue
    fi

    # Subset to ChrX, normalize multiallelics, and index
    bcftools view -S "$ID_FILE" "$VCF_IN" chrX | \
    bcftools norm -m -any -f "$REF" | \
    bcftools annotate --rename-chrs <(echo "chrX chrX") -Oz \
    -o ${OUT_DIR}/TargetALS_chrX_FINAL_${SEX^^}.vcf.gz

    bcftools index -f -t ${OUT_DIR}/TargetALS_chrX_FINAL_${SEX^^}.vcf.gz
done

# --- STEP 3: VERIFICATION ---
echo "----------------------------------------"
echo "Final Sample Counts:"
[ -f ${OUT_DIR}/TargetALS_chrX_FINAL_MALES.vcf.gz ] && \
    echo -n "Males: " && bcftools query -l ${OUT_DIR}/TargetALS_chrX_FINAL_MALES.vcf.gz | wc -l
[ -f ${OUT_DIR}/TargetALS_chrX_FINAL_FEMALES.vcf.gz ] && \
    echo -n "Females: " && bcftools query -l ${OUT_DIR}/TargetALS_chrX_FINAL_FEMALES.vcf.gz | wc -l
echo "----------------------------------------"
