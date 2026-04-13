#!/bin/bash
#SBATCH --job-name=vcf_to_plink_joint
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --time=12:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/vcf_to_plink.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/vcf_to_plink.err

set -euo pipefail

module load PLINK2/avx2_20250707

#----------------------------------------
# PATHS
#----------------------------------------
VCF=/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped.vcf.gz
OUTDIR=/home/zw529/donglab/data/target_ALS/QTL/plink
mkdir -p ${OUTDIR}

RAW_PREFIX=${OUTDIR}/joint_autosomes_raw
QC_SAMPLE_PREFIX=${OUTDIR}/joint_samples_qc
FILTERED_PREFIX=${OUTDIR}/joint_autosomes_filtered
BED_PREFIX=${OUTDIR}/joint_autosomes_filtered_bed
MATRIX_PREFIX=${OUTDIR}/joint_autosomes_matrixEQTL

#----------------------------------------
# STEP 1: VCF → PLINK2 (With Corrected Sex Thresholds)
#----------------------------------------
# Fix: Added specific parameter names for thresholds.
# Step 3 Integration: --chr chr1-22, chrX, chrY ensures we only keep mapped variants.
plink2 --vcf ${VCF} \
       --chr chr1-22, chrX, chrY, chrM \
       --split-par hg38 \
       --impute-sex max-female-xf=0.2 min-male-xf=0.8 \
       --make-pgen \
       --out ${RAW_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# SAMPLE-LEVEL QC (Steps 4, 9, 10)
#----------------------------------------
# Step 4: Check Sex 
# Corrected Syntax: PLINK2 requires the labels max-female-xf and min-male-xf
plink2 --pfile ${RAW_PREFIX} \
       --check-sex max-female-xf=0.2 min-male-xf=0.8 \
       --out ${QC_SAMPLE_PREFIX}

# Step 9: Heterozygosity
plink2 --pfile ${RAW_PREFIX} --het --out ${QC_SAMPLE_PREFIX}

# Step 10: Relatedness
plink2 --pfile ${RAW_PREFIX} --king-cutoff 0.45 --out ${QC_SAMPLE_PREFIX}_relatedness

# --- ID EXTRACTION ---
# Identify Sex Mismatches (Now that thresholds are applied)
grep "PROBLEM" ${QC_SAMPLE_PREFIX}.sexcheck | awk '{print $1, $2}' > ${OUTDIR}/fail_sex.txt || touch ${OUTDIR}/fail_sex.txt

# Identify Heterozygosity Outliers (|F| > 0.2)
awk 'NR>1 && ($6 > 0.2 || $6 < -0.2) {print $1, $2}' ${QC_SAMPLE_PREFIX}.het > ${OUTDIR}/fail_het.txt || touch ${OUTDIR}/fail_het.txt

# Combine all exclusion lists
cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id | sort | uniq > ${OUTDIR}/all_fail_samples.txt

#----------------------------------------
# STEP 2: Filter SNPs & Subjects (Steps 2, 5, 6, 8)
#----------------------------------------
plink2 --pfile ${RAW_PREFIX} \
       --mind 0.05 \
       --geno 0.05 \
       --hwe 1e-6 \
       --maf 0.05 \
       --max-alleles 2 \
       --remove ${OUTDIR}/all_fail_samples.txt \
       --make-pgen \
       --out ${FILTERED_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# OUTPUT GENERATION
#----------------------------------------
plink2 --pfile ${FILTERED_PREFIX} --make-bed --out ${BED_PREFIX}
plink2 --pfile ${FILTERED_PREFIX} --recode A --out ${MATRIX_PREFIX}

echo "QC Pipeline Complete. Final outputs in ${OUTDIR}"
