#!/bin/bash
#SBATCH --job-name=vcf_to_plink_joint
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --time=12:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/vcf_to_plink.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/vcf_to_plink.err

set -euo pipefail

module load PLINK2/avx2_20250707

#----------------------------------------
# INPUT (joint VCF)
#----------------------------------------
VCF=/home/zw529/donglab/data/target_ALS/eQTL/joint_genotyped.vcf.gz

#----------------------------------------
# OUTPUT DIRECTORY
#----------------------------------------
OUTDIR=/home/zw529/donglab/data/target_ALS/eQTL/plink
mkdir -p ${OUTDIR}

#----------------------------------------
# STEP 0: define prefixes
#----------------------------------------
RAW_PREFIX=${OUTDIR}/joint_autosomes_raw
FILTERED_PREFIX=${OUTDIR}/joint_autosomes_filtered
BED_PREFIX=${OUTDIR}/joint_autosomes_filtered_bed
MATRIX_PREFIX=${OUTDIR}/joint_autosomes_matrixEQTL

#----------------------------------------
# STEP 1: VCF → PLINK2 binary (PGEN/PSAM/PVAR) keeping all variants
#----------------------------------------
plink2 --vcf ${VCF} \
       --chr 1-22 \
       --make-pgen \
       --out ${RAW_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# STEP 2: Filter multiallelic + MAF ≥ 0.01 → new PGEN
#----------------------------------------
plink2 --pfile ${RAW_PREFIX} \
       --maf 0.01 \
       --max-alleles 2 \
       --make-pgen \
       --out ${FILTERED_PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# STEP 3: Optional: create BED/BIM/FAM (PLINK1 format)
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

#----------------------------------------
# DONE
#----------------------------------------
echo "Finished PLINK conversion for joint cohort."
echo "All outputs (PGEN, PSAM, PVAR, BED, numeric matrix) are in ${OUTDIR}"
