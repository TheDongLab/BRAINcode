#!/bin/bash
#SBATCH --job-name=TargetALS_Final_Full_Fix
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

mkdir -p ${OUT_DIR}

echo "STARTING FULL WORKFLOW: $(date)"

# ---------------------------------------------------------
# STEP 1: CONVERT VCF TO BINARY (Handle chrX specifically)
# ---------------------------------------------------------
echo "Step 1: Converting VCF to PLINK binary..."
plink --vcf ${VCF_IN} \
      --make-bed \
      --split-x hg38 no-fail \
      --allow-extra-chr \
      --out ${TEMP_PREFIX}

# ---------------------------------------------------------
# STEP 2: SEX AUDIT
# ---------------------------------------------------------
echo "Step 2: Running genetic sex check..."
plink --bfile ${TEMP_PREFIX} \
      --check-sex \
      --allow-extra-chr \
      --out ${OUT_DIR}/final_audit

# ---------------------------------------------------------
# STEP 3: FILTER SAMPLES
# ---------------------------------------------------------
echo "Step 3: Creating clean inclusion lists..."
# SNPSEX 1 = Male, 2 = Female
awk '$4 == 1 {print $1, $2}' ${OUT_DIR}/final_audit.sexcheck > ${OUT_DIR}/ids_male_clean.txt
awk '$4 == 2 {print $1, $2}' ${OUT_DIR}/final_audit.sexcheck > ${OUT_DIR}/ids_female_clean.txt

M_COUNT=$(wc -l < ${OUT_DIR}/ids_male_clean.txt)
F_COUNT=$(wc -l < ${OUT_DIR}/ids_female_clean.txt)
echo "Identified ${M_COUNT} clear males and ${F_COUNT} clear females."

# ---------------------------------------------------------
# STEP 4: EXPORT FINAL CLEAN VCFS (Added --allow-extra-chr)
# ---------------------------------------------------------
echo "Step 4: Exporting final VCFs..."

# Export Clean Males
plink --bfile ${TEMP_PREFIX} \
      --keep ${OUT_DIR}/ids_male_clean.txt \
      --chr X \
      --recode vcf bgz \
      --allow-extra-chr \
      --out ${OUT_DIR}/TargetALS_chrX_CLEAN_MALES

# Export Clean Females
plink --bfile ${TEMP_PREFIX} \
      --keep ${OUT_DIR}/ids_female_clean.txt \
      --chr X \
      --recode vcf bgz \
      --allow-extra-chr \
      --out ${OUT_DIR}/TargetALS_chrX_CLEAN_FEMALES

# ---------------------------------------------------------
# STEP 5: VERIFICATION & CLEANUP
# ---------------------------------------------------------
echo "WORKFLOW COMPLETE: $(date)"
if [[ -f "${OUT_DIR}/TargetALS_chrX_CLEAN_MALES.vcf.gz" && -f "${OUT_DIR}/TargetALS_chrX_CLEAN_FEMALES.vcf.gz" ]]; then
    echo "Success! VCFs generated."
    ls -lh ${OUT_DIR}/TargetALS_chrX_CLEAN_*.vcf.gz
    # rm ${TEMP_PREFIX}.* # Uncomment once you've verified the VCFs
else
    echo "ERROR: VCF generation failed."
    exit 1
fi
