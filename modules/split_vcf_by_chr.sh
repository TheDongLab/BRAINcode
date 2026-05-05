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

# 1. Prep lists
awk 'NR>1 {print $1}' $PSAM_FILE > ${OUTPUT_DIR}/pass_ids.txt
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/full_meta.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="male" {print $1}' > ${OUTPUT_DIR}/males.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="female" {print $1}' > ${OUTPUT_DIR}/females.txt

# ============ APPROACH 1: Standard (Baseline) ============
echo "Running Approach 1..."
APPROACH1_DIR="${OUTPUT_DIR}/approach1_standard"
mkdir -p ${APPROACH1_DIR}

for chr in {1..22} X Y; do
    bcftools view -r chr${chr} -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz
    # Special ploidy fix for chrX only
    if [ "$chr" == "X" ]; then
        python3 << 'EOF'
import gzip, os
vcf = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach1_standard/target_ALS_chrX.vcf.gz"
male_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/males.txt"
PAR = [(10001, 2781479), (155701383, 156030895)]
with open(male_file, 'r') as f: males = set(line.strip() for line in f)
with gzip.open(vcf, 'rt') as inf, gzip.open(vcf+".tmp.gz", 'wt') as outf:
    for line in inf:
        if line.startswith('#'):
            if line.startswith('#CHROM'): samples = line.strip().split('\t')[9:]
            outf.write(line); continue
        cols = line.strip().split('\t'); pos = int(cols[1])
        is_par = any(s <= pos <= e for s, e in PAR)
        for i in range(9, len(cols)):
            if samples[i-9] in males:
                v = cols[i]; a = v[0]; r = v[v.find(":"):] if ":" in v else ""
                cols[i] = f"{a}/{a}{r}" if is_par and "/" not in v[:3] else f"{a}{r}" if not is_par else v
        outf.write('\t'.join(cols) + '\n')
os.replace(vcf+".tmp.gz", vcf)
EOF
    fi
    bcftools index -f -t ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz
done

# ============ APPROACH 2: PAR Separated ============
echo "Running Approach 2 (Linking Autosomes)..."
APPROACH2_DIR="${OUTPUT_DIR}/approach2_par_removed"
mkdir -p ${APPROACH2_DIR}

# Symlink autosomes from App1 to save space
for chr in {1..22} Y; do
    ln -sf ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz ${APPROACH2_DIR}/
    ln -sf ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz.tbi ${APPROACH2_DIR}/
done

# Create the specific X subsets
bcftools view -r chrX:10001-2781479 -S ${OUTPUT_DIR}/pass_ids.txt $INPUT_VCF -Oz -o ${APPROACH2_DIR}/target_ALS_PAR1.vcf.gz
bcftools view -r chrX:155701383-156030895 -S ${OUTPUT_DIR}/pass_ids.txt $INPUT_VCF -Oz -o ${APPROACH2_DIR}/target_ALS_PAR2.vcf.gz
echo -e "chrX\t10001\t2781479\nchrX\t155701383\t156030895" > ${APPROACH2_DIR}/par.bed
bcftools view -T ^${APPROACH2_DIR}/par.bed -r chrX -S ${OUTPUT_DIR}/pass_ids.txt $INPUT_VCF -Oz -o ${APPROACH2_DIR}/target_ALS_chrX_nonPAR.vcf.gz

# Force haploid for males in nonPAR only
python3 << 'EOF'
import gzip, os
vcf = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach2_par_removed/target_ALS_chrX_nonPAR.vcf.gz"
male_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/males.txt"
with open(male_file, 'r') as f: males = set(line.strip() for line in f)
with gzip.open(vcf, 'rt') as inf, gzip.open(vcf+".tmp.gz", 'wt') as outf:
    for line in inf:
        if line.startswith('#'):
            if line.startswith('#CHROM'): samples = line.strip().split('\t')[9:]
            outf.write(line); continue
        cols = line.strip().split('\t')
        for i in range(9, len(cols)):
            if samples[i-9] in males:
                v = cols[i]; a = v[0]; r = v[v.find(":"):] if ":" in v else ""
                cols[i] = f"{a}{r}"
        outf.write('\t'.join(cols) + '\n')
os.replace(vcf+".tmp.gz", vcf)
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
    # Autosomes + X/Y
    for chr in {1..22} X Y; do
        bcftools view -r chr${chr} -S ${SEX_LIST} $INPUT_VCF -Oz -o ${APP3_SEX_DIR}/target_ALS_chr${chr}.vcf.gz
        if [ "$sex" == "males" ] && [ "$chr" == "X" ]; then
            # Apply haploid fix to non-PAR for males
            python3 << 'EOF'
import gzip, os
vcf = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach3_sexstratified/males/target_ALS_chrX.vcf.gz"
PAR = [(10001, 2781479), (155701383, 156030895)]
with gzip.open(vcf, 'rt') as inf, gzip.open(vcf+".tmp.gz", 'wt') as outf:
    for line in inf:
        if line.startswith('#'): outf.write(line); continue
        cols = line.strip().split('\t'); pos = int(cols[1])
        if not any(s <= pos <= e for s, e in PAR):
            for i in range(9, len(cols)):
                v = cols[i]; a = v[0]; r = v[v.find(":"):] if ":" in v else ""
                cols[i] = f"{a}{r}"
        outf.write('\t'.join(cols) + '\n')
os.replace(vcf+".tmp.gz", vcf)
EOF
        fi
        bcftools index -f -t ${APP3_SEX_DIR}/target_ALS_chr${chr}.vcf.gz
    done
done

# ============ VALIDATION SUITE ============
echo -e "\n--- STARTING REFINED VALIDATION ---"

# 1. Check Sample Sex Integrity in Approach 3
for sex in males females; do
    VCF_FILE="${OUTPUT_DIR}/approach3_sexstratified/${sex}/target_ALS_chrX.vcf.gz"
    EXPECTED_IDS="${OUTPUT_DIR}/${sex}.txt"
    ACTUAL_COUNT=$(bcftools query -l $VCF_FILE | wc -l)
    EXPECTED_COUNT=$(wc -l < $EXPECTED_IDS)
    
    echo "Check: Approach 3 ${sex} folder"
    echo "  - Expected: $EXPECTED_COUNT samples | Found: $ACTUAL_COUNT samples"
done

# 2. Check Naming Convention & Header Integrity (Approach 2)
echo -e "\nChecking Chromosome Naming Conventions (Approach 2):"
FILES_TO_CHECK=(
    "${APPROACH2_DIR}/target_ALS_chrX_nonPAR.vcf.gz"
    "${APPROACH2_DIR}/target_ALS_PAR1.vcf.gz"
    "${APPROACH2_DIR}/target_ALS_PAR2.vcf.gz"
    "${OUTPUT_DIR}/approach3_sexstratified/males/target_ALS_chrY.vcf.gz"
)

for f in "${FILES_TO_CHECK[@]}"; do
    if [ -f "$f" ]; then
        BN=$(basename $f)
        CONTIG=$(bcftools view -h "$f" | grep "##contig" | head -n 1)
        COL_CHROM=$(bcftools query -f '%CHROM\n' "$f" | head -n 1)
        echo "File: $BN"
        echo "  - Header Contig: $CONTIG"
        echo "  - Data CHROM Col: $COL_CHROM"
    else
        echo "File: $(basename $f) NOT FOUND!"
    fi
done

# 3. Check Male nonPAR Ploidy (Strictly Haploid Check)
echo -e "\nVerifying nonPAR is ONLY haploid for males..."
# This looks for any genotype containing '/' or '|' in the male-specific or non-PAR files
# We test the Approach 2 nonPAR file (which contains males)
DIPLOID_IN_NONPAR=$(bcftools query -f '[%GT\t]\n' ${APPROACH2_DIR}/target_ALS_chrX_nonPAR.vcf.gz | grep -E "/|\|" | wc -l)

if [ "$DIPLOID_IN_NONPAR" -eq 0 ]; then
    echo "SUCCESS: target_ALS_chrX_nonPAR is 100% haploid for all samples."
else
    echo "ERROR: Found $DIPLOID_IN_NONPAR diploid genotypes in nonPAR file!"
fi

# 4. Confirm PAR regions remain Diploid (They should NOT be haploid)
echo -e "\nConfirming PAR regions remain Diploid..."
PAR1_DIPLOID=$(bcftools query -f '[%GT\t]\n' ${APPROACH2_DIR}/target_ALS_PAR1.vcf.gz | grep -E "/|\|" | wc -l)
if [ "$PAR1_DIPLOID" -gt 0 ]; then
    echo "SUCCESS: PAR1 contains diploid calls (correct for PAR)."
else
    echo "WARNING: PAR1 appears entirely haploid. Check if ploidy fix was over-applied."
fi

echo -e "\n--- ALL VALIDATIONS COMPLETE ---"
