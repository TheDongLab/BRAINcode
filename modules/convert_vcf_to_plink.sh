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
# STEP 0: THE GQ FILTER (Isolated Environment)
#----------------------------------------
# Setting GQ to 0 effectively bypasses the quality filter, 
# allowing all called genotypes to pass into PLINK.
echo "Running Step 1 (GQ Filter 0) in isolated BCFtools environment..."
module --force purge
module load StdEnv
module load GCC/13.3.0
module load BCFtools/1.21-GCC-13.3.0

bcftools filter -e 'FORMAT/GQ < 0' -S . ${VCF_IN} -O z -o ${VCF_FILTERED}
tabix -f -p vcf ${VCF_FILTERED}

#----------------------------------------
# STEP 1: VCF → PLINK2 IMPORT
#----------------------------------------
echo "Switching to PLINK2 environment for remaining QC steps..."
module --force purge
module load StdEnv
module load PLINK2/avx2_20250707
module load R

echo "Importing VCF to PLINK2..."
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

# Cleanup the temp VCF to save space
rm ${VCF_FILTERED} ${VCF_FILTERED}.tbi

#----------------------------------------
# STEP 2: DIAGNOSTICS & HISTOGRAMS
#----------------------------------------
echo "Generating missingness reports..."
plink2 --pfile ${RAW_PREFIX} --missing --out ${OUTDIR}/qc_distribution

Rscript -e "
    # Load data with robust column handling for PLINK2 headers
    smiss <- read.table('${OUTDIR}/qc_distribution.smiss', header=TRUE, check.names=FALSE);
    vmiss <- read.table('${OUTDIR}/qc_distribution.vmiss', header=TRUE, check.names=FALSE);
    
    # Extract values safely
    s_fmiss <- as.numeric(smiss[['F_MISS']]);
    v_fmiss <- as.numeric(vmiss[['F_MISS']]);
    
    pdf('${OUTDIR}/qc_histograms.pdf', width=10, height=5);
    par(mfrow=c(1,2));
    
    # Step 2: Subject Call Rate
    hist(1 - s_fmiss, main='Subject Call Rates (Step 2)', xlab='Call Rate', col='skyblue', breaks=50, xlim=c(0,1));
    abline(v=0.95, col='red', lty=2);
    
    # Step 5: SNP Genotyping Rate
    hist(1 - v_fmiss, main='SNP Genotyping Rates (Step 5)', xlab='Genotyping Rate', col='salmon', breaks=50, xlim=c(0,1));
    abline(v=0.95, col='red', lty=2);
    
    dev.off();
"

#----------------------------------------
# STEP 3: SAMPLE-LEVEL QC
#----------------------------------------
echo "Performing sample-level checks (Sex, Heterozygosity, Relatedness)..."
plink2 --pfile ${RAW_PREFIX} --check-sex max-female-xf=0.2 min-male-xf=0.8 --out ${QC_SAMPLE_PREFIX}
plink2 --pfile ${RAW_PREFIX} --het --out ${QC_SAMPLE_PREFIX}
plink2 --pfile ${RAW_PREFIX} --king-cutoff 0.45 --out ${QC_SAMPLE_PREFIX}_relatedness

# Extract failed IDs
# Use touch/grep logic to ensure scripts don't crash if no failures are found
grep "PROBLEM" ${QC_SAMPLE_PREFIX}.sexcheck | awk '{print $1, $2}' > ${OUTDIR}/fail_sex.txt || touch ${OUTDIR}/fail_sex.txt
awk 'NR>1 && ($6 > 0.2 || $6 < -0.2) {print $1, $2}' ${QC_SAMPLE_PREFIX}.het > ${OUTDIR}/fail_het.txt || touch ${OUTDIR}/fail_het.txt

# Merge all sample failures into one list
cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id | sort | uniq > ${OUTDIR}/all_fail_samples.txt

#----------------------------------------
# STEP 4: FINAL FILTERING
#----------------------------------------
echo "Applying final thresholds from protocol..."
# Note: We apply --geno 0.05 here. Since GQ was 0, this will filter 
# variants based solely on original VCF missingness.
plink2 --pfile ${RAW_PREFIX} \
       --remove ${OUTDIR}/all_fail_samples.txt \
       --mind 0.05 \
       --geno 0.05 \
       --hwe 1e-6 \
       --maf 0.05 \
       --max-alleles 2 \
       --make-pgen \
       --out ${FILTERED_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}
       
#----------------------------------------
# STEP 5 PREP: PRINCIPAL COMPONENT ANALYSIS (PCA)
#----------------------------------------
echo "Running PCA for population stratification..."
plink2 --pfile ${FILTERED_PREFIX} \
       --pca 10 \
       --out ${OUTDIR}/joint_pca

#----------------------------------------
# OUTPUT GENERATION
#----------------------------------------
# 1. Generate Legacy BED files
plink2 --pfile ${FILTERED_PREFIX} --make-bed --out ${BED_PREFIX}

# 2. Generate .raw text matrix (Additive 0/1/2 format) for MatrixEQTL
plink2 --pfile ${FILTERED_PREFIX} --recode A --out ${MATRIX_PREFIX}

echo "QC Pipeline Complete. Output located in ${OUTDIR}"
