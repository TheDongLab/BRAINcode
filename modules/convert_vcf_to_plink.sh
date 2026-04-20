#!/bin/bash
#SBATCH --job-name=vcf_to_plink_joint
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/vcf_to_plink.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/vcf_to_plink.err

set -euo pipefail

#----------------------------------------
# PATHS
#----------------------------------------
VCF_IN=/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz
OUTDIR=/home/zw529/donglab/data/target_ALS/QTL/plink
mkdir -p ${OUTDIR}

VCF_FILTERED=${OUTDIR}/joint_GQ0_temp.vcf.gz
RAW_PREFIX=${OUTDIR}/joint_autosomes_raw
QC_SAMPLE_PREFIX=${OUTDIR}/joint_samples_qc
FILTERED_PREFIX=${OUTDIR}/joint_autosomes_filtered
BED_PREFIX=${OUTDIR}/joint_autosomes_filtered_bed
MATRIX_PREFIX=${OUTDIR}/joint_autosomes_matrixEQTL

#----------------------------------------
# STEP 0: THE GQ FILTER (GQ=0)
#----------------------------------------
echo "Running Step 1 (GQ Filter 0)..."
module --force purge
module load StdEnv
module load GCC/13.3.0
module load BCFtools/1.21-GCC-13.3.0

# Using GQ < 0 to ensure absolutely nothing is filtered at this stage
bcftools filter -e 'FORMAT/GQ < 0' -S . ${VCF_IN} -O z -o ${VCF_FILTERED}
tabix -f -p vcf ${VCF_FILTERED}

#----------------------------------------
# STEP 1: VCF → PLINK2 IMPORT
#----------------------------------------
echo "Importing VCF to PLINK2..."
module --force purge
module load StdEnv
module load PLINK2/avx2_20250707

plink2 --vcf ${VCF_FILTERED} \
       --chr chr1-22, chrX, chrY, chrM \
       --split-par hg38 \
       --set-all-var-ids @:#:\$r:\$a \
       --new-id-max-allele-len 800 \
       --vcf-half-call m \
       --impute-sex max-female-xf=0.2 min-male-xf=0.8 \
       --make-pgen \
       --out ${RAW_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

rm ${VCF_FILTERED} ${VCF_FILTERED}.tbi

#----------------------------------------
# STEP 2: SAMPLE-LEVEL QC
#----------------------------------------
echo "Performing sample-level checks..."
# Check Sex
plink2 --pfile ${RAW_PREFIX} --check-sex max-female-xf=0.2 min-male-xf=0.8 --out ${QC_SAMPLE_PREFIX}
# Check Heterozygosity
plink2 --pfile ${RAW_PREFIX} --het --out ${QC_SAMPLE_PREFIX}
# Check Relatedness
plink2 --pfile ${RAW_PREFIX} --king-cutoff 0.45 --out ${QC_SAMPLE_PREFIX}_relatedness

# Extract failed IDs
# We use '|| true' and 'touch' to prevent the script from exiting if 0 samples fail
grep "PROBLEM" ${QC_SAMPLE_PREFIX}.sexcheck | awk '{print $1, $2}' > ${OUTDIR}/fail_sex.txt || touch ${OUTDIR}/fail_sex.txt
awk 'NR>1 && ($6 > 0.2 || $6 < -0.2) {print $1, $2}' ${QC_SAMPLE_PREFIX}.het > ${OUTDIR}/fail_het.txt || touch ${OUTDIR}/fail_het.txt

# Merge failure lists
# Ensure relatedness output exists before catting
if [ -f ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id ]; then
    cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id | sort | uniq > ${OUTDIR}/all_fail_samples.txt
else
    cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt | sort | uniq > ${OUTDIR}/all_fail_samples.txt
fi

#----------------------------------------
# STEP 3: FINAL QC FILTERING
#----------------------------------------
# 1. Filter out the "sparse" variants first (keep variants present in 95% of people)
plink2 --pfile ${RAW_PREFIX} \
       --geno 0.05 \
       --make-pgen \
       --out ${OUTDIR}/temp_variant_filtered

# 2. NOW check the samples. 
# Because the sparse noise is gone, the samples will suddenly look 99% complete.
plink2 --pfile ${OUTDIR}/temp_variant_filtered \
       --remove ${OUTDIR}/all_fail_samples.txt \
       --mind 0.05 \
       --maf 0.05 \
       --hwe 1e-6 \
       --max-alleles 2 \
       --make-pgen \
       --out ${FILTERED_PREFIX}
       
#----------------------------------------
# STEP 4: PCA & OUTPUT
#----------------------------------------
echo "Running PCA and generating final formats..."
# Top 10 PCs for ancestry covariates
plink2 --pfile ${FILTERED_PREFIX} --pca 10 --out ${OUTDIR}/joint_pca

# Legacy BED and MatrixEQTL .raw
plink2 --pfile ${FILTERED_PREFIX} --make-bed --out ${BED_PREFIX}
plink2 --pfile ${FILTERED_PREFIX} --recode A --nonfounders --out ${MATRIX_PREFIX}   # ensures we get hard-call counts (0, 1, 2) and NOT probabilities (0.0 - 1.0)

echo "Pipeline Finished."
