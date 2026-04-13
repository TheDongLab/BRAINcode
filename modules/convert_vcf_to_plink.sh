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
module load R

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
# STEP 1: FORCE VCF → PLINK2 IMPORT
#----------------------------------------
# --allow-no-sex is the key here. It stops the crash.
# --split-par hg38 handles the PAR regions correctly for hg38 coordinates.
plink2 --vcf ${VCF} \
       --chr chr1-22, chrX, chrY, chrM \
       --split-par hg38 \
       --allow-no-sex \
       --make-pgen \
       --out ${RAW_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# STEP 2 & 5: THE "FEEL" REPORTS (Histograms)
#----------------------------------------
# We do this immediately so you can see why the 95% filter is failing.
plink2 --pfile ${RAW_PREFIX} --missing --out ${OUTDIR}/qc_distribution

Rscript -e "
    smiss <- read.table('${OUTDIR}/qc_distribution.smiss', header=TRUE);
    vmiss <- read.table('${OUTDIR}/qc_distribution.vmiss', header=TRUE);
    pdf('${OUTDIR}/qc_histograms.pdf', width=10, height=5);
    par(mfrow=c(1,2));
    hist(1 - smiss\$F_MISS, main='Subject Call Rates', xlab='Call Rate', col='skyblue', breaks=50);
    abline(v=0.95, col='red', lty=2);
    hist(1 - vmiss\$F_MISS, main='SNP Genotyping Rates', xlab='Genotyping Rate', col='salmon', breaks=50);
    abline(v=0.95, col='red', lty=2);
    dev.off();
"

#----------------------------------------
# STEP 4: IMPUTE & CHECK SEX (The logic you needed)
#----------------------------------------
# 1. First, we impute sex and create a NEW .psam file.
plink2 --pfile ${RAW_PREFIX} \
       --impute-sex max-female-xf=0.2 min-male-xf=0.8 \
       --make-psam \
       --out ${OUTDIR}/imputed_metadata

# 2. IMPORTANT: We replace the "empty" psam with the "imputed" one.
# This "teaches" PLINK the sex of your samples.
cp ${OUTDIR}/imputed_metadata.psam ${RAW_PREFIX}.psam

# 3. Now run the check-sex and other checks.
plink2 --pfile ${RAW_PREFIX} --check-sex max-female-xf=0.2 min-male-xf=0.8 --out ${QC_SAMPLE_PREFIX}
plink2 --pfile ${RAW_PREFIX} --het --out ${QC_SAMPLE_PREFIX}
plink2 --pfile ${RAW_PREFIX} --king-cutoff 0.45 --out ${QC_SAMPLE_PREFIX}_relatedness

#----------------------------------------
# EXTRACTION & FINAL FILTER
#----------------------------------------
grep "PROBLEM" ${QC_SAMPLE_PREFIX}.sexcheck | awk '{print \$1, \$2}' > ${OUTDIR}/fail_sex.txt || touch ${OUTDIR}/fail_sex.txt
awk 'NR>1 && (\$6 > 0.2 || \$6 < -0.2) {print \$1, \$2}' ${QC_SAMPLE_PREFIX}.het > ${OUTDIR}/fail_het.txt || touch ${OUTDIR}/fail_het.txt

cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id | sort | uniq > ${OUTDIR}/all_fail_samples.txt

# Final Step: If this still results in "No samples remaining", your histograms 
# from the R script will show you exactly how much to lower the 0.05 thresholds.
plink2 --pfile ${RAW_PREFIX} \
       --remove ${OUTDIR}/all_fail_samples.txt \
       --mind 0.05 \
       --geno 0.05 \
       --hwe 1e-6 \
       --maf 0.05 \
       --make-pgen \
       --out ${FILTERED_PREFIX}

#----------------------------------------
# OUTPUT GENERATION
#----------------------------------------
plink2 --pfile ${FILTERED_PREFIX} --make-bed --out ${BED_PREFIX}
plink2 --pfile ${FILTERED_PREFIX} --recode A --out ${MATRIX_PREFIX}

echo "QC Pipeline Complete. Final outputs in ${OUTDIR}"
