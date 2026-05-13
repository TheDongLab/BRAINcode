#!/bin/bash
#SBATCH --job-name=TargetALS_Full_Audit
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/audit_full_%j.out

# Keep running even if individual commands fail
set -uo pipefail

VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
PSAM="/home/zw529/donglab/data/target_ALS/QTL/diagnostics/perfect_match.psam"
FAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/diagnostics/perfect_match.fam"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"

# ---------------------------------------------------------
# STEP 1: PLINK 1.9 - The "Legacy" Tests
# ---------------------------------------------------------
echo "Running PLINK 1.9 Comprehensive Audit..."
module load PLINK/1.9b_7.11-x86_64  || echo "P1.9 Load Failed"

# We chain these with '|| true' so one failure doesn't kill the batch
plink --vcf ${VCF_IN} --make-bed --split-x hg38 no-fail --out ${OUT_DIR}/p19_temp || true

# Apply HH fix and run full diagnostics
plink --bfile ${OUT_DIR}/p19_temp --chr X --set-hh-missing --make-bed --out ${OUT_DIR}/p19_chrX_clean || true

# Full Suite: Sex Check, Missingness, Freq, Heterozygosity, and HWE
plink --bfile ${OUT_DIR}/p19_chrX_clean \
      --check-sex \
      --missing \
      --freq \
      --het \
      --hardy \
      --out ${OUT_DIR}/p19_full_audit || echo "P1.9 Diagnostics encountered an error"

# ---------------------------------------------------------
# STEP 2: PLINK 2.0 - The "Modern" & High-Speed Audit
# ---------------------------------------------------------
echo "Running PLINK 2.0 Comprehensive Audit..."
module load PLINK2/avx2_20250707 || echo "P2.0 Load Failed"

# PLINK 2 diagnostics are much faster and handle VCF metadata better
plink2 --vcf ${VCF_IN} \
       --psam ${PSAM} \
       --split-par hg38 \
       --check-sex \
       --missing \
       --freq \
       --het \
       --hardy \
       --sample-diff \
       --out ${OUT_DIR}/p20_full_audit || echo "P2.0 Diagnostics encountered an error"

# ---------------------------------------------------------
# STEP 3: Cleanup
# ---------------------------------------------------------
rm -f ${OUT_DIR}/p19_temp.*
echo "All audits complete. Check ${OUT_DIR} for .sexcheck, .het, .hardy, and .imiss files."
