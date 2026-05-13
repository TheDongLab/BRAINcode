#!/bin/bash
#SBATCH --job-name=TargetALS_Final_Full_Pipeline
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/FULL_PIPELINE_%j.out

set -uo pipefail

# --- CONFIGURATION ---
VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"
TEMP_PREFIX="${OUT_DIR}/tmp_p19_process"

# Load modules
module load PLINK/1.90-beta6.10
module load bcftools/1.17

mkdir -p ${OUT_DIR}

echo "STARTING FULL WORKFLOW: $(date)"

# ---------------------------------------------------------
# STEP 1: CONVERT VCF TO BINARY (Handle chrX specifically)
# ---------------------------------------------------------
echo "Step 1: Converting VCF to PLINK binary and splitting PAR..."
plink --vcf ${VCF_IN} \
      --make-bed \
      --split-x hg38 no-fail \
      --allow-extra-chr \
      --out ${TEMP_PREFIX}

# ---------------------------------------------------------
# STEP 2: SEX AUDIT (Determine who is who genetically)
# ---------------------------------------------------------
echo "Step 2: Running genetic sex check..."
plink --bfile ${TEMP_PREFIX} \
      --check-sex \
      --allow-extra-chr \
      --out ${OUT_DIR}/final_audit

# ---------------------------------------------------------
# STEP 3: FILTER SAMPLES (Remove the 68 Ambiguous Samples)
# ---------------------------------------------------------
echo "Step 3: Creating clean exclusion/inclusion lists..."

# Extract IDs based on SNPSEX (Column 4 of .sexcheck)
# SNPSEX 1 = Male, 2 = Female, 0 = Ambiguous
awk '$4 == 1 {print $1, $2}' ${OUT_DIR}/final_audit.sexcheck > ${OUT_DIR}/ids_male_clean.txt
awk '$4 == 2 {print $1, $2}' ${OUT_DIR}/final_audit.sexcheck > ${OUT_DIR}/ids_female_clean.txt

# Count the results for the log
M_COUNT=$(wc -l < ${OUT_DIR}/ids_male_clean.txt)
F_COUNT=$(wc -l < ${OUT_DIR}/ids_female_clean.txt)
echo "Found ${M_COUNT} clear males and ${F_COUNT} clear females."

# ---------------------------------------------------------
# STEP 4: EXPORT FINAL CLEAN VCFS (Stratified by Sex)
# ---------------------------------------------------------
echo "Step 4: Exporting final VCFs for TOPMed upload..."

# Clean Males ChrX
plink --bfile ${TEMP_PREFIX} \
      --keep ${OUT_DIR}/ids_male_clean.txt \
      --chr X \
      --recode vcf bgz \
      --out ${OUT_DIR}/TargetALS_chrX_CLEAN_MALES

# Clean Females ChrX
plink --bfile ${TEMP_PREFIX} \
      --keep ${OUT_DIR}/ids_female_clean.txt \
      --chr X \
      --recode vcf bgz \
      --out ${OUT_DIR}/TargetALS_chrX_CLEAN_FEMALES

# ---------------------------------------------------------
# STEP 5: CLEAN UP INTERMEDIATE FILES
# ---------------------------------------------------------
echo "Step 5: Cleaning up temp binary files..."
rm ${TEMP_PREFIX}.bed ${TEMP_PREFIX}.bim ${TEMP_PREFIX}.fam ${TEMP_PREFIX}.log ${TEMP_PREFIX}.nosex

echo "WORKFLOW COMPLETE: $(date)"
echo "Files ready for upload in ${OUT_DIR}:"
ls -lh ${OUT_DIR}/TargetALS_chrX_CLEAN_*.vcf.gz
