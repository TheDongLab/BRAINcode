#!/bin/bash
#SBATCH --job-name=split_vcf_by_chr
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

export LC_ALL=C
module load BCFtools/1.21

# Paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
BASE_OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/chromosome_joint_vcfs"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"

mkdir -p ${OUTPUT_DIR}

# 1. Prep lists - REFINED TO RESCUE XX/XY FROM COLUMN 28
awk 'NR>1 {print $1}' ${PSAM_FILE} > ${OUTPUT_DIR}/pass_ids.txt

tr -d '\r' < ${METADATA} | awk -F',' 'NR>1 {
    s=tolower($0); 
    # Check standard sex column string OR column 28 for XX/XY
    if(s ~ /female/ || $28 == "XX") sex="female"; 
    else if(s ~ /male/ || $28 == "XY") sex="male"; 
    else sex="unknown"; 
    print $2, sex
}' | sort -u > ${OUTPUT_DIR}/full_meta.txt

grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="male" {print $1}' > ${OUTPUT_DIR}/males.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="female" {print $1}' > ${OUTPUT_DIR}/females.txt

# ============ APPROACH 1: Standard (Baseline) ============
echo "Running Approach 1..."
APPROACH1_DIR="${OUTPUT_DIR}/approach1_standard"
mkdir -p ${APPROACH1_DIR}

for chr in {1..22} X Y; do
    bcftools view -r chr${chr} -S ${OUTPUT_DIR}/pass_ids.txt --force-samples ${INPUT_VCF} -Oz -o ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz
    if [ "$chr" == "X" ]; then
        export VCF_PROC="${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz"
        export MALE_PROC="${OUTPUT_DIR}/males.txt"
        python3 << 'EOF'
import gzip, os, subprocess
vcf = os.environ['VCF_PROC']
male_file = os.environ['MALE_PROC']
PAR = [(10001, 2781479), (155701383, 156030895)]
with open(male_file, 'r') as f: males = set(line.strip() for line in f)

with gzip.open(vcf, 'rt') as inf, open(vcf + ".tmp", 'w') as outf:
    for line in inf:
        if line.startswith('#'):
            if line.startswith('#CHROM'): samples = line.strip().split('\t')[9:]
            outf.write(line); continue
        cols = line.strip().split('\t'); pos = int(cols[1])
        is_par = any(s <= pos <= e for s, e in PAR)
        for i in range(9, len(cols)):
            v = cols[i]
            if v.startswith("."): continue
            parts = v.split(':'); gt = parts[0]; rest = ":" + ":".join(parts[1:]) if len(parts) > 1 else ""
            a = gt[0]
            if samples[i-9] in males:
                new_gt = f"{a}/{a}" if is_par and "/" not in gt and "|" not in gt else a if not is_par else gt
            else:
                new_gt = f"{a}/{a}" if "/" not in gt and "|" not in gt else gt
            cols[i] = f"{new_gt}{rest}"
        outf.write('\t'.join(cols) + '\n')
subprocess.run(["bcftools", "view", vcf + ".tmp", "-Oz", "-o", vcf], check=True)
os.remove(vcf + ".tmp")
EOF
    fi
    bcftools index -f -t ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz
done

# ============ APPROACH 2: PAR Separated ============
echo "Running Approach 2..."
APPROACH2_DIR="${OUTPUT_DIR}/approach2_par_removed"
mkdir -p ${APPROACH2_DIR}

for chr in {1..22} Y; do
    ln -sf ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz ${APPROACH2_DIR}/
    ln -sf ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz.tbi ${APPROACH2_DIR}/
done

bcftools view -r chrX:10001-2781479 -S ${OUTPUT_DIR}/pass_ids.txt ${INPUT_VCF} -Oz -o ${APPROACH2_DIR}/target_ALS_PAR1.vcf.gz
bcftools view -r chrX:155701383-156030895 -S ${OUTPUT_DIR}/pass_ids.txt ${INPUT_VCF} -Oz -o ${APPROACH2_DIR}/target_ALS_PAR2.vcf.gz
echo -e "chrX\t10001\t2781479\nchrX\t155701383\t156030895" > ${APPROACH2_DIR}/par.bed
bcftools view -T ^${APPROACH2_DIR}/par.bed -r chrX -S ${OUTPUT_DIR}/pass_ids.txt ${INPUT_VCF} -Oz -o ${APPROACH2_DIR}/target_ALS_chrX_nonPAR.vcf.gz

export MALE_PROC="${OUTPUT_DIR}/males.txt"
export VCF_PROC="${APPROACH2_DIR}/target_ALS_chrX_nonPAR.vcf.gz"
python3 << 'EOF'
import gzip, os, subprocess
vcf = os.environ['VCF_PROC']
male_file = os.environ['MALE_PROC']
with open(male_file, 'r') as f: males = set(line.strip() for line in f)

with gzip.open(vcf, 'rt') as inf, open(vcf + ".tmp", 'w') as outf:
    for line in inf:
        if line.startswith('#'):
            if line.startswith('#CHROM'): samples = line.strip().split('\t')[9:]
            outf.write(line); continue
        cols = line.strip().split('\t')
        for i in range(9, len(cols)):
            v = cols[i]
            if v.startswith("."): continue
            parts = v.split(':'); gt = parts[0]; rest = ":" + ":".join(parts[1:]) if len(parts) > 1 else ""
            a = gt[0]
            new_gt = a if samples[i-9] in males else (f"{a}/{a}" if "/" not in gt and "|" not in gt else gt)
            cols[i] = f"{new_gt}{rest}"
        outf.write('\t'.join(cols) + '\n')
subprocess.run(["bcftools", "view", vcf + ".tmp", "-Oz", "-o", vcf], check=True)
os.remove(vcf + ".tmp")
EOF
bcftools index -f -t ${APPROACH2_DIR}/target_ALS_chrX_nonPAR.vcf.gz
bcftools index -f -t ${APPROACH2_DIR}/target_ALS_PAR1.vcf.gz
bcftools index -f -t ${APPROACH2_DIR}/target_ALS_PAR2.vcf.gz

# ============ APPROACH 3: Sex-stratified (Strict Biallelic) ============
echo "Running Approach 3..."
for sex in males females; do
    APP3_SEX_DIR="${OUTPUT_DIR}/approach3_sexstratified_fixed/${sex}"
    mkdir -p ${APP3_SEX_DIR}
    SEX_LIST="${OUTPUT_DIR}/${sex}.txt"
    
    for chr in {1..22} X Y; do
        OUT_VCF="${APP3_SEX_DIR}/target_ALS_chr${chr}.vcf.gz"
        echo "Processing Chromosome ${chr} for ${sex}..."
        
        # KEY CHANGE: Remove -m2 -M2 -v snps filter — it kills haploids
        bcftools view -r chr${chr} -S ${SEX_LIST} "${INPUT_VCF}" -Ou | \
        bcftools norm -m-any -Ou | \
        bcftools norm -d both -Oz -o "${OUT_VCF}"

        if [ "$chr" == "X" ] && [ "$sex" == "males" ]; then
            export CURRENT_VCF="${OUT_VCF}"
            python3 << 'EOF'
import gzip, os, subprocess
vcf = os.environ['CURRENT_VCF']
PAR = [(10001, 2781479), (155701383, 156030895)]

with open(vcf + ".tmp", 'w') as outf, gzip.open(vcf, 'rt') as inf:
    for line in inf:
        if line.startswith('#'):
            outf.write(line); continue
        cols = line.strip().split('\t')
        pos = int(cols[1])
        is_par = any(s <= pos <= e for s, e in PAR)
        for i in range(9, len(cols)):
            v = cols[i]
            if v.startswith("."): continue
            parts = v.split(':')
            gt = parts[0]
            suffix = ":" + ":".join(parts[1:]) if len(parts) > 1 else ""
            a = gt[0]
            new_gt = f"{a}/{a}{suffix}" if is_par else f"{a}{suffix}"
            cols[i] = new_gt
        outf.write('\t'.join(cols) + '\n')
subprocess.run(["bcftools", "view", vcf + ".tmp", "-Oz", "-o", vcf], check=True)
os.remove(vcf + ".tmp")
EOF
        fi
        bcftools index -f -t "${OUT_VCF}"
    done
done

# ============ VALIDATION SUITE ============
echo -e "\n--- STARTING FINAL VALIDATION ---"
for sex in males females; do
    VCF_FILE="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach3_sexstratified/${sex}/target_ALS_chrX.vcf.gz"
    ACTUAL_COUNT=$(bcftools query -l ${VCF_FILE} | wc -l)
    GHOSTS=$(bcftools query -f '[%GT\t]\n' ${VCF_FILE} | grep -E "[2-9]" | wc -l)
    echo "Check: Approach 3 ${sex} - Samples: ${ACTUAL_COUNT}, Ghost Alleles: ${GHOSTS}"
done

FEMALE_HAP_IN_X=$(bcftools query -f '[%GT\t]\n' /home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach3_sexstratified/females/target_ALS_chrX.vcf.gz | tr '\t' '\n' | grep -v "\." | grep -E "^[01]$" | wc -l)
if [ "${FEMALE_HAP_IN_X}" -eq 0 ]; then echo "SUCCESS: Female X 100% Diploid"; else echo "FAIL: ${FEMALE_HAP_IN_X} haploids in Female X"; fi
