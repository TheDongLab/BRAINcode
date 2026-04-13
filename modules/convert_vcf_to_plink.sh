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
# STEP 1: VCF → PLINK2 (MANDATORY IMPORT LOGIC)
#----------------------------------------
# PLINK2 requires --impute-sex HERE to handle chrX during import.
# Step 3: Satisfied by limiting to chr1-22, X, Y, M (primary hg38 assembly).
plink2 --vcf ${VCF} \
       --chr chr1-22, chrX, chrY, chrM \
       --split-par hg38 \
       --impute-sex max-female-xf=0.2 min-male-xf=0.8 \
       --make-pgen \
       --out ${RAW_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# STEP 2 & 5: THE "FEEL" REPORTS (Histograms)
#----------------------------------------
plink2 --pfile ${RAW_PREFIX} --missing --out ${OUTDIR}/qc_distribution

Rscript -e "
    smiss <- read.table('${OUTDIR}/qc_distribution.smiss', header=TRUE);
    vmiss <- read.table('${OUTDIR}/qc_distribution.vmiss', header=TRUE);
    pdf('${OUTDIR}/qc_histograms.pdf', width=10, height=5);
    par(mfrow=c(1,2));
    # Step 2: Subject Call Rate
    hist(1 - smiss\$F_MISS, main='Subject Call Rates', xlab='Call Rate', col='skyblue', breaks=50);
    abline(v=0.95, col='red', lty=2);
    # Step 5: SNP Genotyping Rate
    hist(1 - vmiss\$F_MISS, main='SNP Genotyping Rates', xlab='Genotyping Rate', col='salmon', breaks=50);
    abline(v=0.95, col='red', lty=2);
    dev.off();
"

#----------------------------------------
# STEP 4, 9, 10: QC CHECKS
#----------------------------------------
# Step 4: Check Sex (Thresholds: 0.2/0.8)
plink2 --pfile ${RAW_PREFIX} \
       --check-sex max-female-xf=0.2 min-male-xf=0.8 \
       --out ${QC_SAMPLE_PREFIX}

# Step 9: Heterozygosity
plink2 --pfile ${RAW_PREFIX} --het --out ${QC_SAMPLE_PREFIX}

# Step 10: Relatedness (KING PI_HAT > 0.9)
plink2 --pfile ${RAW_PREFIX} --king-cutoff 0.45 --out ${QC_SAMPLE_PREFIX}_relatedness

# --- ID EXTRACTION ---
grep "PROBLEM" ${QC_SAMPLE_PREFIX}.sexcheck | awk '{print $1, $2}' > ${OUTDIR}/fail_sex.txt || touch ${OUTDIR}/fail_sex.txt
awk 'NR>1 && ($6 > 0.2 || $6 < -0.2) {print $1, $2}' ${QC_SAMPLE_PREFIX}.het > ${OUTDIR}/fail_het.txt || touch ${OUTDIR}/fail_het.txt
cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id | sort | uniq > ${OUTDIR}/all_fail_samples.txt

#----------------------------------------
# STEP 2, 5, 6, 8: FINAL FILTERING
#----------------------------------------
# Step 6: HWE (catch artifacts/mishaps)
# Step 8: MAF 0.05
plink2 --pfile ${RAW_PREFIX} \
       --remove ${OUTDIR}/all_fail_samples.txt \
       --mind 0.05 \
       --geno 0.05 \
       --hwe 1e-6 \
       --maf 0.05 \
       --max-alleles 2 \
       --make-pgen \
       --out ${FILTERED_PREFIX}

#----------------------------------------
# OUTPUT GENERATION
#----------------------------------------
plink2 --pfile ${FILTERED_PREFIX} --make-bed --out ${BED_PREFIX}
plink2 --pfile ${FILTERED_PREFIX} --recode A --out ${MATRIX_PREFIX}

echo "QC Pipeline Complete. Final outputs in ${OUTDIR}"
