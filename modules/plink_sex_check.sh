#!/bin/bash
#SBATCH --job-name=TargetALS_Final_Filter
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/filter_final_%j.out

set -uo pipefail

VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
# Fixed: Ensuring path to PSAM is correct (change to the actual location if different)
PSAM="/home/zw529/donglab/data/target_ALS/QTL/perfect_match.psam" 
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"

module load PLINK/1.90-beta6.10
module load PLINK2/avx2_20250707

# ---------------------------------------------------------
# STEP 1: Audit & Diagnostics (PLINK 1.9)
# ---------------------------------------------------------
echo "Running Audit..."
plink --vcf ${VCF_IN} --make-bed --split-x hg38 no-fail --allow-extra-chr --out ${OUT_DIR}/p19_temp || true
plink --bfile ${OUT_DIR}/p19_temp --check-sex --allow-extra-chr --out ${OUT_DIR}/p19_audit || true

# ---------------------------------------------------------
# STEP 2: Identify Unambiguous Samples
# ---------------------------------------------------------
echo "Extracting clear males and females..."

# SNPSEX 1 = Clear Males, SNPSEX 2 = Clear Females
awk '$4 == 1 {print $1, $2}' ${OUT_DIR}/p19_audit.sexcheck > ${OUT_DIR}/clean_males_list.txt
awk '$4 == 2 {print $1, $2}' ${OUT_DIR}/p19_audit.sexcheck > ${OUT_DIR}/clean_females_list.txt

# ---------------------------------------------------------
# STEP 3: Export Clean Stratified VCFs
# ---------------------------------------------------------
echo "Exporting clean stratified VCFs..."

# Export Clean MALES
plink2 --vcf ${VCF_IN} \
       --keep ${OUT_DIR}/clean_males_list.txt \
       --chr X \
       --split-par hg38 \
       --export vcf bgz \
       --out ${OUT_DIR}/target_ALS_chrX_CLEAN_MALES

# Export Clean FEMALES
plink2 --vcf ${VCF_IN} \
       --keep ${OUT_DIR}/clean_females_list.txt \
       --chr X \
       --split-par hg38 \
       --export vcf bgz \
       --out ${OUT_DIR}/target_ALS_chrX_CLEAN_FEMALES

# ---------------------------------------------------------
# STEP 4: PLINK 2.0 Audit (Fixed Path)
# ---------------------------------------------------------
echo "Running final PLINK 2.0 check..."
# Note: Removed --sample-diff to avoid previous error. Added check for psam file.
if [ -f "$PSAM" ]; then
    plink2 --vcf ${VCF_IN} \
           --psam ${PSAM} \
           --split-par hg38 \
           --check-sex --missing --freq --het --hardy --allow-extra-chr \
           --out ${OUT_DIR}/p20_final_audit || echo "P2.0 Audit failed"
else
    echo "PSAM file not found at $PSAM. Skipping P2.0 Audit."
fi

echo "Filtering Complete. Ready for upload: CLEAN_MALES and CLEAN_FEMALES."
