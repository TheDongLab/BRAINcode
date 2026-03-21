#!/bin/bash
#SBATCH --job-name=vcf_to_plink_joint
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
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

PREFIX=${OUTDIR}/joint_autosomes

#----------------------------------------
# STEP 1: VCF → PLINK BED (autosomes only, biallelic, MAF filter)
#----------------------------------------
plink2 --vcf ${VCF} \
       --chr 1-22 \  # ignore sex chrs
       --maf 0.01 \
       --max-alleles 2 \
       --make-bed \
       --out ${PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# STEP 2: PLINK → numeric allele matrix for MatrixEQTL
#----------------------------------------
plink2 --bfile ${PREFIX} \
       --recode A \
       --out ${PREFIX}.matrixEQTL \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# DONE
#----------------------------------------
echo "Finished PLINK conversion (autosomes only) for joint cohort"
