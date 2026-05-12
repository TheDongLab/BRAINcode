#!/bin/bash
#SBATCH --job-name=split_vcf_by_chr
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

export LC_ALL=C
module load BCFtools/1.21

INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
BASE_OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"

# Define the subject to exclude due to sex mismatch
EXCLUDE_ID="NEUTN269WMP"

mkdir -p ${BASE_OUTPUT_DIR}/chromosome_joint_vcfs
OUTPUT_DIR="${BASE_OUTPUT_DIR}/chromosome_joint_vcfs"

# Generate pass_ids and explicitly filter out the mismatch subject
awk 'NR>1 {print $1}' $PSAM_FILE | grep -v "$EXCLUDE_ID" > ${OUTPUT_DIR}/pass_ids.txt

# Create metadata mappings
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/full_meta.txt

# These lists will now be clean because they reference the filtered pass_ids.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="male" {print $1}' > ${OUTPUT_DIR}/males.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="female" {print $1}' > ${OUTPUT_DIR}/females.txt

echo "QC samples: $(wc -l < ${OUTPUT_DIR}/pass_ids.txt), Males: $(wc -l < ${OUTPUT_DIR}/males.txt), Females: $(wc -l < ${OUTPUT_DIR}/females.txt)"

# ============ APPROACH 1: Standard ============
echo "APPROACH 1: Standard separation"
APPROACH1_DIR="${OUTPUT_DIR}/approach1_standard"
mkdir -p ${APPROACH1_DIR}
bcftools view -r chrX -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${APPROACH1_DIR}/target_ALS_chrX.vcf.gz

python3 << 'EOF'
import gzip, os
vcf_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach1_standard/target_ALS_chrX.vcf.gz"
male_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/males.txt"
temp_file = vcf_file + ".tmp.gz"
PAR = [(10001, 2781479), (155701383, 156030895)]
with open(male_file, 'r') as f:
    males = set(line.strip() for line in f)
with gzip.open(vcf_file, 'rt') as inf, gzip.open(temp_file, 'wt') as outf:
    for line in inf:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
            outf.write(line)
            continue
        cols = line.strip().split('\t')
        pos = int(cols[1])
        is_par = any(start <= pos <= end for start, end in PAR)
        for i in range(9, len(cols)):
            if samples[i-9] in males:
                val = cols[i]
                allele = val[0]
                sep_idx = val.find(":")
                rest = val[sep_idx:] if sep_idx != -1 else ""
                if is_par:
                    if "/" not in val[:3] and "|" not in val[:3]:
                        cols[i] = f"{allele}/{allele}{rest}"
                else:
                    cols[i] = f"{allele}{rest}"
        outf.write('\t'.join(cols) + '\n')
os.replace(temp_file, vcf_file)
EOF

bcftools view ${APPROACH1_DIR}/target_ALS_chrX.vcf.gz -Oz -o ${APPROACH1_DIR}/target_ALS_chrX_final.vcf.gz
mv ${APPROACH1_DIR}/target_ALS_chrX_final.vcf.gz ${APPROACH1_DIR}/target_ALS_chrX.vcf.gz
bcftools index -f -t ${APPROACH1_DIR}/target_ALS_chrX.vcf.gz

for chr in {1..22} Y; do
    bcftools view -r chr${chr} -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${APPROACH1_DIR}/target_ALS_chr${chr}.vcf.gz
done

# ============ APPROACH 2: PAR Removed ============
echo "APPROACH 2: PAR as separate chromosomes"
APPROACH2_DIR="${OUTPUT_DIR}/approach2_par_removed"
mkdir -p ${APPROACH2_DIR}
bcftools view -r chrX -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${APPROACH2_DIR}/target_ALS_chrX_full.vcf.gz

python3 << 'PYEOF'
import gzip, os
vcf_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach2_par_removed/target_ALS_chrX_full.vcf.gz"
male_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/males.txt"
nonpar_vcf = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach2_par_removed/target_ALS_chrX_nonPAR.vcf.gz"
par1_vcf = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach2_par_removed/target_ALS_PAR1.vcf.gz"
par2_vcf = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach2_par_removed/target_ALS_PAR2.vcf.gz"
PAR1 = (10001, 2781479)
PAR2 = (155701383, 156030895)
with open(male_file, 'r') as f:
    males = set(line.strip() for line in f)
header_lines = []
sample_indices = {}
with gzip.open(vcf_file, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            header_lines.append(line)
            samples = line.strip().split('\t')[9:]
            sample_indices = {s: i for i, s in enumerate(samples)}
            break
        else:
            header_lines.append(line)
with gzip.open(vcf_file, 'rt') as inf:
    with gzip.open(nonpar_vcf, 'wt') as nonpar_f, gzip.open(par1_vcf, 'wt') as par1_f, gzip.open(par2_vcf, 'wt') as par2_f:
        for header_line in header_lines:
            nonpar_f.write(header_line)
            par1_f.write(header_line)
            par2_f.write(header_line)
        for line in inf:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            pos = int(cols[1])
            in_par1 = PAR1[0] <= pos <= PAR1[1]
            in_par2 = PAR2[0] <= pos <= PAR2[1]
            in_nonpar = not (in_par1 or in_par2)
            if in_nonpar:
                for i in range(9, len(cols)):
                    if samples[i-9] in males:
                        val = cols[i]
                        allele = val[0]
                        sep_idx = val.find(":")
                        rest = val[sep_idx:] if sep_idx != -1 else ""
                        cols[i] = f"{allele}{rest}"
                nonpar_f.write('\t'.join(cols) + '\n')
            if in_par1:
                cols[0] = 'PAR1'
                par1_f.write('\t'.join(cols) + '\n')
            if in_par2:
                cols[0] = 'PAR2'
                par2_f.write('\t'.join(cols) + '\n')
os.remove(vcf_file)
PYEOF

for vcf in ${APPROACH2_DIR}/target_ALS_chrX_nonPAR.vcf.gz ${APPROACH2_DIR}/target_ALS_PAR1.vcf.gz ${APPROACH2_DIR}/target_ALS_PAR2.vcf.gz; do
    [ -f "$vcf" ] && bcftools view ${vcf} -Oz -o ${vcf}.tmp && mv ${vcf}.tmp ${vcf} && bcftools index -f -t ${vcf}
done

bcftools view -r chrY -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${APPROACH2_DIR}/target_ALS_chrY_full.vcf.gz

python3 << 'PYEOF'
import gzip, os
vcf_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach2_par_removed/target_ALS_chrY_full.vcf.gz"
male_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/males.txt"
nonpar_vcf = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/approach2_par_removed/target_ALS_chrY_nonPAR.vcf.gz"
PAR1 = (10001, 2781479)
PAR2 = (155701383, 156030895)
with open(male_file, 'r') as f:
    males = set(line.strip() for line in f)
header_lines = []
with gzip.open(vcf_file, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            header_lines.append(line)
            samples = line.strip().split('\t')[9:]
            break
        else:
            header_lines.append(line)
with gzip.open(vcf_file, 'rt') as inf:
    with gzip.open(nonpar_vcf, 'wt') as nonpar_f:
        for header_line in header_lines:
            nonpar_f.write(header_line)
        for line in inf:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            pos = int(cols[1])
            in_par = (PAR1[0] <= pos <= PAR1[1]) or (PAR2[0] <= pos <= PAR2[1])
            if not in_par:
                for i in range(9, len(cols)):
                    if samples[i-9] in males:
                        val = cols[i]
                        allele = val[0]
                        sep_idx = val.find(":")
                        rest = val[sep_idx:] if sep_idx != -1 else ""
                        cols[i] = f"{allele}{rest}"
                nonpar_f.write('\t'.join(cols) + '\n')
os.remove(vcf_file)
PYEOF

bcftools view ${APPROACH2_DIR}/target_ALS_chrY_nonPAR.vcf.gz -Oz -o ${APPROACH2_DIR}/target_ALS_chrY_nonPAR_final.vcf.gz
mv ${APPROACH2_DIR}/target_ALS_chrY_nonPAR_final.vcf.gz ${APPROACH2_DIR}/target_ALS_chrY_nonPAR.vcf.gz
bcftools index -f -t ${APPROACH2_DIR}/target_ALS_chrY_nonPAR.vcf.gz

for chr in {1..22}; do
    bcftools view -r chr${chr} -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${APPROACH2_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${APPROACH2_DIR}/target_ALS_chr${chr}.vcf.gz
done

# ============ APPROACH 3: Sex-stratified ============
echo "APPROACH 3: Sex-stratified (males/females separate)"
for sex in males females; do
    SEX_LIST_FILE="${OUTPUT_DIR}/${sex}.txt"
    APPROACH3_SEX_DIR="${OUTPUT_DIR}/approach3_sexstratified/${sex}"
    mkdir -p ${APPROACH3_SEX_DIR}
    bcftools view -r chrX -S ${SEX_LIST_FILE} --force-samples $INPUT_VCF -Oz -o ${APPROACH3_SEX_DIR}/target_ALS_chrX.vcf.gz
    python3 << EOF
import gzip, os
vcf_file = "${APPROACH3_SEX_DIR}/target_ALS_chrX.vcf.gz"
is_male = ('${sex}' == 'males')
temp_file = vcf_file + ".tmp.gz"
PAR = [(10001, 2781479), (155701383, 156030895)]
with gzip.open(vcf_file, 'rt') as inf, gzip.open(temp_file, 'wt') as outf:
    for line in inf:
        if line.startswith('#'):
            outf.write(line)
            continue
        cols = line.strip().split('\t')
        pos = int(cols[1])
        is_par = any(start <= pos <= end for start, end in PAR)
        if is_male and not is_par:
            for i in range(9, len(cols)):
                val = cols[i]
                allele = val[0]
                sep_idx = val.find(":")
                rest = val[sep_idx:] if sep_idx != -1 else ""
                cols[i] = f"{allele}{rest}"
        outf.write('\t'.join(cols) + '\n')
os.replace(temp_file, vcf_file)
EOF
    bcftools view ${APPROACH3_SEX_DIR}/target_ALS_chrX.vcf.gz -Oz -o ${APPROACH3_SEX_DIR}/target_ALS_chrX_final.vcf.gz
    mv ${APPROACH3_SEX_DIR}/target_ALS_chrX_final.vcf.gz ${APPROACH3_SEX_DIR}/target_ALS_chrX.vcf.gz
    bcftools index -f -t ${APPROACH3_SEX_DIR}/target_ALS_chrX.vcf.gz
    for chr in {1..22} Y; do
        bcftools view -r chr${chr} -S ${SEX_LIST_FILE} --force-samples $INPUT_VCF -Oz -o ${APPROACH3_SEX_DIR}/target_ALS_chr${chr}.vcf.gz
        bcftools index -f -t ${APPROACH3_SEX_DIR}/target_ALS_chr${chr}.vcf.gz
    done
done

echo "All approaches complete: approach1_standard, approach2_par_removed, approach3_sexstratified"
