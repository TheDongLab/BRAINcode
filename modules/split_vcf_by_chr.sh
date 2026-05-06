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

# 1. Prep lists - Refined to catch 'male/female' regardless of column shifts
awk 'NR>1 {print $1}' $PSAM_FILE > ${OUTPUT_DIR}/pass_ids.txt
tr -d '\r' < $METADATA | awk -F',' 'NR>1 {s=tolower($0); if(s ~ /female/) sex="female"; else if(s ~ /male/) sex="male"; else sex="unknown"; print $2, sex}' | sort -u > ${OUTPUT_DIR}/full_meta.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="male" {print $1}' > ${OUTPUT_DIR}/males.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="female" {print $1}' > ${OUTPUT_DIR}/females.txt

# ============ APPROACH 1: Standard (Baseline) ============
echo "Running Approach 1..."
APPROACH1_DIR="${OUTPUT_DIR}/approach1_standard"
mkdir -p ${APPROACH1_DIR}

for chr in {1..22} X Y; do
    bcftools view -r chr${chr} -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz
    if [ "$chr" == "X" ]; then
        python3 << 'EOF'
import gzip, os, subprocess
vcf = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach1_standard/target_ALS_chrX.vcf.gz"
male_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/males.txt"
PAR = [(10001, 2781479), (155701383, 156030895)]
with open(male_file, 'r') as f: males = set(line.strip() for line in f)
with gzip.open(vcf, 'rt') as inf, open(vcf+".tmp", 'w') as outf:
    for line in inf:
        if line.startswith('#'):
            if line.startswith('#CHROM'): samples = line.strip().split('\t')[9:]
            outf.write(line); continue
        cols = line.strip().split('\t'); pos = int(cols[1])
        is_par = any(s <= pos <= e for s, e in PAR)
        for i in range(9, len(cols)):
            v = cols[i]
            if v.startswith("."): continue
            a = v[0]; r = v[v.find(":"):] if ":" in v else ""
            if samples[i-9] in males:
                cols[i] = f"{a}/{a}{r}" if is_par and "/" not in v[:3] else f"{a}{r}" if not is_par else v
            else: # Force Females to Diploid
                if "/" not in v[:3] and "|" not in v[:3]: cols[i] = f"{a}/{a}{r}"
        outf.write('\t'.join(cols) + '\n')
# Re-compress with bcftools to ensure BGZF format for indexing
subprocess.run(["bcftools", "view", vcf+".tmp", "-Oz", "-o", vcf], check=True)
os.remove(vcf+".tmp")
EOF
    fi
    bcftools index -f -t ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz
done

# ============ APPROACH 2: PAR Separated ============
echo "Running Approach 2 (Linking Autosomes)..."
APPROACH2_DIR="${OUTPUT_DIR}/approach2_par_removed"
mkdir -p ${APPROACH2_DIR}

for chr in {1..22} Y; do
    ln -sf ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz ${APPROACH2_DIR}/
    ln -sf ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz.tbi ${APPROACH2_DIR}/
done

bcftools view -r chrX:10001-2781479 -S ${OUTPUT_DIR}/pass_ids.txt $INPUT_VCF -Oz -o ${APPROACH2_DIR}/target_ALS_PAR1.vcf.gz
bcftools view -r chrX:155701383-156030895 -S ${OUTPUT_DIR}/pass_ids.txt $INPUT_VCF -Oz -o ${APPROACH2_DIR}/target_ALS_PAR2.vcf.gz
echo -e "chrX\t10001\t2781479\nchrX\t155701383\t156030895" > ${APPROACH2_DIR}/par.bed
bcftools view -T ^${APPROACH2_DIR}/par.bed -r chrX -S ${OUTPUT_DIR}/pass_ids.txt $INPUT_VCF -Oz -o ${APPROACH2_DIR}/target_ALS_chrX_nonPAR.vcf.gz

export MALE_FILE_APP2="${OUTPUT_DIR}/males.txt"
export VCF_APP2="${APPROACH2_DIR}/target_ALS_chrX_nonPAR.vcf.gz"
python3 << 'EOF'
import gzip, os, subprocess
vcf = os.environ['VCF_APP2']
male_file = os.environ['MALE_FILE_APP2']
with open(male_file, 'r') as f: males = set(line.strip() for line in f)
with gzip.open(vcf, 'rt') as inf, open(vcf+".tmp", 'w') as outf:
    for line in inf:
        if line.startswith('#'):
            if line.startswith('#CHROM'): samples = line.strip().split('\t')[9:]
            outf.write(line); continue
        cols = line.strip().split('\t')
        for i in range(9, len(cols)):
            v = cols[i]
            if v.startswith("."): continue
            a = v[0]; r = v[v.find(":"):] if ":" in v else ""
            if samples[i-9] in males:
                cols[i] = f"{a}{r}"
            else: # Force Females to Diploid
                if "/" not in v[:3] and "|" not in v[:3]: cols[i] = f"{a}/{a}{r}"
        outf.write('\t'.join(cols) + '\n')
subprocess.run(["bcftools", "view", vcf+".tmp", "-Oz", "-o", vcf], check=True)
os.remove(vcf+".tmp")
EOF
bcftools index -f -t ${APPROACH2_DIR}/target_ALS_chrX_nonPAR.vcf.gz
bcftools index -f -t ${APPROACH2_DIR}/target_ALS_PAR1.vcf.gz
bcftools index -f -t ${APPROACH2_DIR}/target_ALS_PAR2.vcf.gz

# ============ APPROACH 3: Sex-stratified ============
echo "Running Approach 3..."
for sex in males females; do
    APP3_SEX_DIR="${OUTPUT_DIR}/approach3_sexstratified/${sex}"
    mkdir -p ${APP3_SEX_DIR}
    SEX_LIST="${OUTPUT_DIR}/${sex}.txt"
    
    for chr in {1..22} X Y; do
        # We add 'bcftools norm -m-any' to split multi-allelic sites 
        # and 'view -m2 -M2' to ensure only biallelic SNPs proceed.
        bcftools view -r chr${chr} -S ${SEX_LIST} $INPUT_VCF | \
        bcftools norm -m-any --fasta-ref $REF_FASTA | \
        bcftools view -m2 -M2 -v snps -Oz -o ${APP3_SEX_DIR}/target_ALS_chr${chr}.vcf.gz
        
        if [ "$chr" == "X" ]; then
            export CURRENT_VCF="${APP3_SEX_DIR}/target_ALS_chr${chr}.vcf.gz"
            export CURRENT_SEX="$sex"
            python3 << 'EOF'
import gzip, os, subprocess

vcf = os.environ['CURRENT_VCF']
sex = os.environ['CURRENT_SEX']
# GRCh38 PAR regions for X
PAR = [(10001, 2781479), (155701383, 156030895)]

with gzip.open(vcf, 'rt') as inf, open(vcf+".tmp", 'w') as outf:
    for line in inf:
        if line.startswith('#'):
            outf.write(line)
            continue
        
        cols = line.strip().split('\t')
        pos = int(cols[1])
        is_par = any(s <= pos <= e for s, e in PAR)
        
        for i in range(9, len(cols)):
            v = cols[i]
            if v.startswith("."): 
                continue
            
            # Split the genotype from other fields (AD, DP, etc.)
            parts = v.split(':')
            gt = parts[0]
            rest = ":" + ":".join(parts[1:]) if len(parts) > 1 else ""
            
            # Extract first allele safely
            a1 = gt[0]
            
            if sex == "males":
                if not is_par:
                    # Force haploid (0 or 1)
                    new_gt = a1
                else:
                    # Force diploid in PAR (0/0, 0/1, 1/1)
                    if "/" not in gt and "|" not in gt:
                        new_gt = f"{a1}/{a1}"
                    else:
                        new_gt = gt
            else:
                # Force Females to Diploid everywhere
                if "/" not in gt and "|" not in gt:
                    new_gt = f"{a1}/{a1}"
                else:
                    new_gt = gt
            
            cols[i] = f"{new_gt}{rest}"
            
        outf.write('\t'.join(cols) + '\n')

subprocess.run(["bcftools", "view", vcf+".tmp", "-Oz", "-o", vcf], check=True)
os.remove(vcf+".tmp")
EOF
        fi
        bcftools index -f -t ${APP3_SEX_DIR}/target_ALS_chr${chr}.vcf.gz
    done
done

# ============ VALIDATION SUITE ============
echo -e "\n--- STARTING REFINED VALIDATION ---"

for sex in males females; do
    VCF_FILE="${OUTPUT_DIR}/approach3_sexstratified/${sex}/target_ALS_chrX.vcf.gz"
    ACTUAL_COUNT=$(bcftools query -l $VCF_FILE | wc -l)
    echo "Check: Approach 3 ${sex} folder - Found: $ACTUAL_COUNT samples"
    
    # NEW: Check for ghost/multi-allelic alleles (anything > 1)
    GHOST_ALLELES=$(bcftools query -f '[%GT\t]\n' $VCF_FILE | grep -E "[2-9]" | wc -l)
    if [ "$GHOST_ALLELES" -eq 0 ]; then
        echo "  - Biallelic Check: SUCCESS (No alleles > 1)"
    else
        echo "  - Biallelic Check: ERROR (Found $GHOST_ALLELES genotypes with alleles > 1!)"
    fi
done

echo -e "\nVerifying Female ChrX is ONLY diploid/missing..."
# Improved check to ensure no single-digit genotypes (0 or 1) exist
FEMALE_HAP_IN_X=$(bcftools query -f '[%GT\t]\n' ${OUTPUT_DIR}/approach3_sexstratified/females/target_ALS_chrX.vcf.gz | \
    tr '\t' '\n' | grep -v "\." | grep -E "^[01]$" | wc -l)

if [ "$FEMALE_HAP_IN_X" -eq 0 ]; then
    echo "SUCCESS: Female ChrX is 100% diploid/missing."
else
    echo "ERROR: Found $FEMALE_HAP_IN_X haploid calls in female file!"
fi

echo -e "\n--- ALL VALIDATIONS COMPLETE ---"
