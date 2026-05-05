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
# PATHS
#----------------------------------------
# NEW: Path for the filtered imputed VCF
VCF_IN=/home/zw529/donglab/data/target_ALS/QTL/target_ALS_imputed_filtered_joint.vcf.gz

# OBSOLETE/NORMAL VCF PATH:
# VCF_IN=/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz

OUTDIR=/home/zw529/donglab/data/target_ALS/QTL/plink
mkdir -p ${OUTDIR}

# These are now only needed if running Step 0 (Normal VCF)
VCF_FILTERED=${OUTDIR}/joint_GQ0_temp.vcf.gz
VCF_BIALLELIC=${OUTDIR}/joint_biallelic_temp.vcf.gz

RAW_PREFIX=${OUTDIR}/joint_autosomes_raw
QC_SAMPLE_PREFIX=${OUTDIR}/joint_samples_qc
FILTERED_PREFIX=${OUTDIR}/joint_autosomes_filtered
BED_PREFIX=${OUTDIR}/joint_autosomes_filtered_bed
MATRIX_PREFIX=${OUTDIR}/joint_autosomes_matrixEQTL

#----------------------------------------
# STEP 0: PRE-PROCESSING (OBSOLETE FOR IMPUTED)
#----------------------------------------
# Note: Imputed data from TOPMed is already normalized and filtered.
# Only uncomment these steps if you are switching back to "Normal" VCF.

# echo "Running Step 0a (GQ Filter)..."
# module --force purge
# module load StdEnv
# module load GCC/13.3.0
# module load BCFtools/1.21-GCC-13.3.0
# bcftools filter -e 'FORMAT/GQ < 0' -S . ${VCF_IN} -O z -o ${VCF_FILTERED}
# tabix -f -p vcf ${VCF_FILTERED}

# echo "Running Step 0b (Multiallelic to Biallelic)..."
# bcftools norm -m -any ${VCF_FILTERED} -O z -o ${VCF_BIALLELIC}
# tabix -f -p vcf ${VCF_BIALLELIC}
# rm ${VCF_FILTERED} ${VCF_FILTERED}.tbi

#----------------------------------------
# STEP 1: VCF → PLINK2 IMPORT
#----------------------------------------
echo "Importing VCF to PLINK2..."
module --force purge
module load StdEnv
module load PLINK2/avx2_20250707

# NEW: Import directly from the Imputed VCF_IN
plink2 --vcf ${VCF_IN} \
       --chr chr1-22, chrX, chrY, chrM \
       --split-par hg38 \
       --set-all-var-ids @:#:\$r:\$a \
       --new-id-max-allele-len 800 \
       --vcf-half-call m \
       --impute-sex max-female-xf=0.2 min-male-xf=0.8 \
       --make-pgen \
       --out ${RAW_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

# OBSOLETE IMPORT (Used the temp biallelic file):
# plink2 --vcf ${VCF_BIALLELIC} ...
# rm ${VCF_BIALLELIC} ${VCF_BIALLELIC}.tbi

#----------------------------------------
# STEP 2: SAMPLE-LEVEL QC
#----------------------------------------
echo "Performing sample-level checks..."
plink2 --pfile ${RAW_PREFIX} --check-sex max-female-xf=0.2 min-male-xf=0.8 --out ${QC_SAMPLE_PREFIX}
plink2 --pfile ${RAW_PREFIX} --het --out ${QC_SAMPLE_PREFIX}
plink2 --pfile ${RAW_PREFIX} --king-cutoff 0.45 --out ${QC_SAMPLE_PREFIX}_relatedness

grep "PROBLEM" ${QC_SAMPLE_PREFIX}.sexcheck | awk '{print $1, $2}' > ${OUTDIR}/fail_sex.txt || touch ${OUTDIR}/fail_sex.txt
awk 'NR>1 && ($6 > 0.2 || $6 < -0.2) {print $1, $2}' ${QC_SAMPLE_PREFIX}.het > ${OUTDIR}/fail_het.txt || touch ${OUTDIR}/fail_het.txt

if [ -f ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id ]; then
    cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id | sort | uniq > ${OUTDIR}/all_fail_samples.txt
else
    cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt | sort | uniq > ${OUTDIR}/all_fail_samples.txt
fi

#----------------------------------------
# STEP 3: FINAL QC FILTERING
#----------------------------------------
# Note: MAF/R2 filters were already applied to the Imputed VCF, 
# but these steps clean up sample-level dropouts.

plink2 --pfile ${RAW_PREFIX} \
       --geno 0.05 \
       --make-pgen \
       --out ${OUTDIR}/temp_variant_filtered

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
plink2 --pfile ${FILTERED_PREFIX} --pca 10 --out ${OUTDIR}/joint_pca
plink2 --pfile ${FILTERED_PREFIX} --make-bed --out ${BED_PREFIX}

# Force hard-call counts (0, 1, 2)
plink2 --pfile ${FILTERED_PREFIX} \
       --recode A \
       --hard-call-threshold 0.49 \
       --nonfounders \
       --out ${MATRIX_PREFIX}

echo "Pipeline Finished."
