#!/bin/bash
#SBATCH --job-name=vcf_to_plink_joint
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH -p day
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/vcf_to_plink.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/vcf_to_plink.err

set -euo pipefail

module load PLINK2/avx2_20250707
module load Python/3.12.3-GCCcore-13.3.0

#----------------------------------------
# INPUT
#----------------------------------------
VCF=/home/zw529/donglab/data/target_ALS/eQTL/joint_genotyped.vcf.gz

BASE=/home/zw529/donglab/data/target_ALS
EQTL=${BASE}/eQTL
OUTDIR=${EQTL}/plink
mkdir -p ${OUTDIR}

cd ${EQTL}

PREFIX=${OUTDIR}/joint
SEXFILE=${EQTL}/sex.txt
MULTIALLELIC_LIST=${OUTDIR}/joint_multiallelic.txt

#----------------------------------------
# STEP 0: auto-generate sex.txt from metadata
#----------------------------------------
python3 - <<EOF
import pandas as pd

meta = pd.read_csv("${BASE}/targetALS_wgs_metadata.csv")
meta.columns = meta.columns.str.strip()

IID_col = "Externalsampleid"
SEX_col = "Sex"

def map_sex(x):
    x=str(x).strip().lower()
    if x in ["male","m","1"]:
        return "1"
    elif x in ["female","f","2"]:
        return "2"
    else:
        return "NA"

meta["SEX_PLINK"] = meta[SEX_col].apply(map_sex)

out = meta[[IID_col,"SEX_PLINK"]].dropna()
out.to_csv("${SEXFILE}", sep=" ", index=False, header=["IID","SEX"])
EOF

#----------------------------------------
# STEP 1a: list multiallelic variants for documentation
#----------------------------------------
plink2 --vcf ${VCF} \
       --max-alleles 2 list \
       --out ${MULTIALLELIC_LIST}

#----------------------------------------
# STEP 1b: VCF → BED (remove multiallelics, filter MAF, handle sex chromosomes)
#----------------------------------------
plink2 --vcf ${VCF} \
       --max-alleles 2 strict \
       --maf 0.01 \
       --impute-sex \
       --split-par b38 \
       --make-bed \
       --out ${PREFIX} \
       --threads ${SLURM_CPUS_PER_TASK}

#----------------------------------------
# STEP 2: numeric allele matrix for MatrixEQTL
#----------------------------------------
plink2 --bfile ${PREFIX} \
       --recode A \
       --out ${PREFIX}.matrixEQTL \
       --threads ${SLURM_CPUS_PER_TASK}

echo "Finished PLINK conversion for joint cohort"
