#!/bin/bash
#SBATCH --job-name=TargetALS_Unified_X_Pipeline
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/UNIFIED_X_PIPELINE_%j.out

set -uo pipefail

# --- CONFIGURATION ---
VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"
TEMP_PREFIX="${OUT_DIR}/tmp_p19_process"

# Load modules
module load PLINK/1.90-beta6.10
module load BCFtools/1.21

mkdir -p ${OUT_DIR}

echo "STARTING UNIFIED WORKFLOW: $(date)"

# ---------------------------------------------------------
# STEP 1: SUBJECT-LEVEL QC (Filter for the 413 subjects)
# ---------------------------------------------------------
echo "Step 1: Extracting the 413 subjects that passed Subject-QC..."

# Get IDs from your validated .psam file (excluding the first column/header if needed)
# This handles the "14 subjects that failed" logic by relying on your pre-filtered PSAM
awk 'NR > 1 {print $1, $2}' ${PSAM_FILE} > ${OUT_DIR}/subject_qc_pass_list.txt

# ---------------------------------------------------------
# STEP 2: GENETIC AUDIT (On the QC-passed subset)
# ---------------------------------------------------------
echo "Step 2: Converting VCF and running Genetic Sex Audit..."

plink --vcf ${VCF_IN} \
      --keep ${OUT_DIR}/subject_qc_pass_list.txt \
      --make-bed \
      --split-x hg38 no-fail \
      --allow-extra-chr \
      --out ${TEMP_PREFIX}

plink --bfile ${TEMP_PREFIX} \
      --check-sex \
      --allow-extra-chr \
      --out ${OUT_DIR}/genetic_sex_audit

# ---------------------------------------------------------
# STEP 3: INTERSECT FILTERS (Final Inclusion Lists)
# ---------------------------------------------------------
echo "Step 3: Creating final inclusion lists (Unambiguous Genetic Sex + Subject-QC)..."

# Extract IDs where SNPSEX is clear (1 or 2)
awk '$4 == 1 {print $1, $2}' ${OUT_DIR}/genetic_sex_audit.sexcheck > ${OUT_DIR}/final_males.txt
awk '$4 == 2 {print $1, $2}' ${OUT_DIR}/genetic_sex_audit.sexcheck > ${OUT_DIR}/final_females.txt

echo "Final Counts: $(wc -l < ${OUT_DIR}/final_males.txt) Males, $(wc -l < ${OUT_DIR}/final_females.txt) Females."

# ---------------------------------------------------------
# STEP 4: EXPORT & ENCODING FIX (Rename 23 to chrX)
# ---------------------------------------------------------
echo "Step 4: Exporting VCFs and fixing hg38 encoding..."

# Create renaming map for bcftools
echo "23 chrX" > ${OUT_DIR}/rename_chrs.txt

for SEX in males females; do
    UPPER_SEX=${SEX^^}
    
    # Export using PLINK 1.9 (Temporary intermediate)
    plink --bfile ${TEMP_PREFIX} \
          --keep ${OUT_DIR}/final_${SEX}.txt \
          --chr X \
          --recode vcf bgz \
          --allow-extra-chr \
          --out ${OUT_DIR}/tmp_export_${SEX}

    # Use BCFtools to rename '23' to 'chrX' and finalize
    bcftools annotate --rename-chrs ${OUT_DIR}/rename_chrs.txt \
        ${OUT_DIR}/tmp_export_${SEX}.vcf.gz \
        -Oz -o ${OUT_DIR}/TargetALS_chrX_FINAL_${UPPER_SEX}.vcf.gz

    # Index for TOPMed
    bcftools index -t ${OUT_DIR}/TargetALS_chrX_FINAL_${UPPER_SEX}.vcf.gz
done

# ---------------------------------------------------------
# STEP 5: CLEANUP
# ---------------------------------------------------------
echo "Step 5: Cleaning up intermediate files..."
rm ${OUT_DIR}/tmp_export_*
rm ${TEMP_PREFIX}.*
rm ${OUT_DIR}/rename_chrs.txt

echo "WORKFLOW COMPLETE: $(date)"
ls -lh ${OUT_DIR}/TargetALS_chrX_FINAL_*.vcf.gz
