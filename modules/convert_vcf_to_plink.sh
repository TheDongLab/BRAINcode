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
# INPUT (joint VCF)
#----------------------------------------
VCF=/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped.vcf.gz

#----------------------------------------
# OUTPUT DIRECTORY
#----------------------------------------
OUTDIR=/home/zw529/donglab/data/target_ALS/QTL/plink
mkdir -p ${OUTDIR}

#----------------------------------------
# STEP 0: define prefixes
#----------------------------------------
RAW_PREFIX=${OUTDIR}/joint_autosomes_raw
QC_SAMPLE_PREFIX=${OUTDIR}/joint_samples_qc
FILTERED_PREFIX=${OUTDIR}/joint_autosomes_filtered
BED_PREFIX=${OUTDIR}/joint_autosomes_filtered_bed
MATRIX_PREFIX=${OUTDIR}/joint_autosomes_matrixEQTL

#----------------------------------------
# STEP 1: VCF → PLINK2 binary
#----------------------------------------
plink2 --vcf ${VCF} \
       --chr 1-22 \
       --make-pgen \
       --out ${RAW_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# NEW: SAMPLE-LEVEL QC (Steps 4, 9, 10)
#----------------------------------------
# Step 4: Check Sex (Note: Requires X chromosome; if autosomes only, this results in a report)
plink2 --pfile ${RAW_PREFIX} --check-sex --out ${QC_SAMPLE_PREFIX}

# Step 9: Heterozygosity (Inbreeding coefficient F)
plink2 --pfile ${RAW_PREFIX} --het --out ${QC_SAMPLE_PREFIX}

# Step 10: Relatedness (KING-robust) - Filters pairs with PI_HAT > 0.9
# Note: --make-king-table generates the list; --king-cutoff removes one of each pair
plink2 --pfile ${RAW_PREFIX} \
       --king-cutoff 0.45 \
       --out ${QC_SAMPLE_PREFIX}_relatedness

#----------------------------------------
# STEP 2: Filter SNPs & Subjects (Steps 2, 5, 6, 8)
#----------------------------------------
# Integrating:
# Step 2: --mind 0.05 (Subject call rate < 95%) [cite: 5]
# Step 5: --geno 0.05 (SNP genotype rate < 95%) [cite: 15]
# Step 6: --hwe 1e-6 (Hardy-Weinberg Equilibrium) [cite: 16]
# Step 8: --maf 0.05 (Minor Allele Frequency) [cite: 20]
plink2 --pfile ${RAW_PREFIX} \
       --mind 0.05 \
       --geno 0.05 \
       --hwe 1e-6 \
       --maf 0.05 \
       --max-alleles 2 \
       --remove ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id \
       --make-pgen \
       --out ${FILTERED_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# STEP 3: Create BED/BIM/FAM
#----------------------------------------
plink2 --pfile ${FILTERED_PREFIX} \
       --make-bed \
       --out ${BED_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# STEP 4: Numeric allele matrix for MatrixEQTL
#----------------------------------------
plink2 --pfile ${FILTERED_PREFIX} \
       --recode A \
       --out ${MATRIX_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

echo "Finished QC and PLINK conversion."
