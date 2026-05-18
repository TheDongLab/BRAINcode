#!/bin/bash
#SBATCH --job-name=vcf_to_plink_imputed
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/vcf_to_plink.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/vcf_to_plink.err

set -euo pipefail

#----------------------------------------
# PATHS & PREFIXES
#----------------------------------------
VCF_IN=/home/zw529/donglab/data/target_ALS/QTL/target_ALS_imputed_filtered_joint.vcf.gz
OUTDIR=/home/zw529/donglab/data/target_ALS/QTL/plink
mkdir -p ${OUTDIR}

# Temp files for Step 0 (if ever uncommented)
VCF_FILTERED=${OUTDIR}/joint_GQ0_temp.vcf.gz
VCF_BIALLELIC=${OUTDIR}/joint_biallelic_temp.vcf.gz

# --- OLD AUTOSOME-ONLY PREFIXES (HASHED OUT) ---
# RAW_PREFIX=${OUTDIR}/joint_autosomes_raw
# QC_SAMPLE_PREFIX=${OUTDIR}/joint_samples_qc
# FILTERED_PREFIX=${OUTDIR}/joint_autosomes_filtered
# BED_PREFIX=${OUTDIR}/joint_autosomes_filtered_bed
# MATRIX_PREFIX=${OUTDIR}/joint_autosomes_matrixEQTL

# --- NEW ALL-CHROMOSOME PREFIXES ---
RAW_PREFIX=${OUTDIR}/joint_all_chrs_raw
QC_SAMPLE_PREFIX=${OUTDIR}/joint_samples_qc
FILTERED_PREFIX=${OUTDIR}/joint_all_chrs_filtered
BED_PREFIX=${OUTDIR}/joint_all_chrs_filtered_bed
MATRIX_PREFIX=${OUTDIR}/joint_all_chrs_matrixEQTL

#----------------------------------------
# STEP 0: PRE-PROCESSING (OBSOLETE FOR IMPUTED VCFs)
#----------------------------------------
# Only uncomment if switching back to "Normal" raw VCFs.
#
# echo "Running Step 0a (GQ Filter)..."
# module --force purge
# module load StdEnv
# module load BCFtools/1.21
# bcftools filter -e 'FORMAT/GQ < 0' -S . ${VCF_IN} -O z -o ${VCF_FILTERED}
# bcftools index -t ${VCF_FILTERED}
#
# echo "Running Step 0b (Multiallelic to Biallelic)..."
# bcftools norm -m -any ${VCF_FILTERED} -O z -o ${VCF_BIALLELIC}
# bcftools index -t ${VCF_BIALLELIC}
# rm ${VCF_FILTERED} ${VCF_FILTERED}.tbi
# VCF_IN=${VCF_BIALLELIC} 

#----------------------------------------
# STEP 1: VCF → PLINK2 IMPORT
#----------------------------------------
module --force purge
module load StdEnv
module load PLINK2/avx2_20250707

#================================================================================
# VERSION 1: WITH SEX CHROMOSOMES (ACTIVE)
#================================================================================
echo "Importing VCF with Sex Chromosomes..."
plink2 --vcf ${VCF_IN} \
        --chr chr1-22,chrX,chrY,chrM \
        --split-par hg38 \
        --set-all-var-ids @:#:\$r:\$a \
        --new-id-max-allele-len 800 \
        --vcf-half-call m \
        --impute-sex max-female-xf=0.2 min-male-xf=0.8 \
        --make-pgen \
        --out ${RAW_PREFIX} \
        --threads ${SLURM_CPUS_PER_TASK}

echo "Running sex-check validation..."
plink2 --pfile ${RAW_PREFIX} --check-sex max-female-xf=0.2 min-male-xf=0.8 --out ${QC_SAMPLE_PREFIX}
grep "PROBLEM" ${QC_SAMPLE_PREFIX}.sexcheck | awk '{print $1, $2}' > ${OUTDIR}/fail_sex.txt || touch ${OUTDIR}/fail_sex.txt


#================================================================================
# VERSION 2: WITHOUT SEX CHROMOSOMES (HASHED OUT)
#================================================================================
# echo "Importing Autosome-only VCF..."
# plink2 --vcf ${VCF_IN} \
#        --chr chr1-22 \
#        --set-all-var-ids @:#:\$r:\$a \
#        --new-id-max-allele-len 800 \
#        --vcf-half-call m \
#        --make-pgen \
#        --out ${RAW_PREFIX} \
#        --threads ${SLURM_CPUS_PER_TASK}
#
# touch ${OUTDIR}/fail_sex.txt

#----------------------------------------
# SAMPLE-LEVEL QC (Shared Logic)
#----------------------------------------
echo "Performing sample-level checks..."
plink2 --pfile ${RAW_PREFIX} --het --out ${QC_SAMPLE_PREFIX}
plink2 --pfile ${RAW_PREFIX} --king-cutoff 0.45 --out ${QC_SAMPLE_PREFIX}_relatedness

awk 'NR>1 && ($6 > 0.2 || $6 < -0.2) {print $1, $2}' ${QC_SAMPLE_PREFIX}.het > ${OUTDIR}/fail_het.txt || touch ${OUTDIR}/fail_het.txt

if [ -f ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id ]; then
    cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id | sort | uniq > ${OUTDIR}/all_fail_samples.txt
else
    cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt | sort | uniq > ${OUTDIR}/all_fail_samples.txt
fi

#----------------------------------------
# FINAL QC FILTERING
#----------------------------------------
echo "Applying variant missingness filter (--geno 0.05)..."
plink2 --pfile ${RAW_PREFIX} \
       --geno 0.05 \
       --make-pgen \
       --out ${OUTDIR}/temp_variant_filtered

echo "Applying final hard sample filters, MAF, and HWE thresholds..."
plink2 --pfile ${OUTDIR}/temp_variant_filtered \
       --remove ${OUTDIR}/all_fail_samples.txt \
       --mind 0.05 \
       --maf 0.05 \
       --hwe 1e-6 \
       --max-alleles 2 \
       --make-pgen \
       --out ${FILTERED_PREFIX}
       
#----------------------------------------
# PCA & OUTPUT
#----------------------------------------
echo "Running PCA and generating final formats..."
plink2 --pfile ${FILTERED_PREFIX} --pca 10 --out ${OUTDIR}/joint_pca
plink2 --pfile ${FILTERED_PREFIX} --make-bed --out ${BED_PREFIX}

plink2 --pfile ${FILTERED_PREFIX} \
       --recode A \
       --hard-call-threshold 0.49 \
       --nonfounders \
       --out ${MATRIX_PREFIX}

echo "Pipeline Finished."
