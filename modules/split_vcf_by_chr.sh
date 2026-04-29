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
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"

mkdir -p $OUTPUT_DIR

# --- STEP 1: PREP LISTS ---
awk 'NR>1 {print $1}' $PSAM_FILE > ${OUTPUT_DIR}/pass_ids.txt
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/full_meta.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="male" {print $1}' > ${OUTPUT_DIR}/males.txt

# --- STEP 2: GENERATE FILES (STANDARD) ---
# We generate chrX first so we can patch it immediately
echo "Extracting chrX..."
bcftools view -r chrX -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz || echo "Warning: chrX extraction encountered an error"

# Generate autosomes and Y
for chr in {1..22} Y; do
    echo "Extracting chr${chr}..."
    bcftools view -r chr${chr} -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

# --- STEP 3: PATCH CHRX IN-PLACE ---
echo "Applying Python ploidy patch to target_ALS_chrX.vcf.gz..."
python3 << 'EOF'
import gzip
import os

vcf_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/target_ALS_chrX.vcf.gz"
male_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/males.txt"
temp_file = vcf_file + ".tmp.gz"

with open(male_file, 'r') as f:
    males = set(line.strip() for line in f)

try:
    with gzip.open(vcf_file, 'rt') as inf, gzip.open(temp_file, 'wt') as outf:
        for line in inf:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    samples = line.strip().split('\t')[9:]
                outf.write(line)
                continue
            
            cols = line.strip().split('\t')
            for i in range(9, len(cols)):
                if samples[i-9] in males:
                    # Strip the second allele and separator: e.g., '0/0:AD...' -> '0:AD...'
                    val = cols[i]
                    colon_idx = val.find(':')
                    if colon_idx != -1:
                        cols[i] = val[0] + val[colon_idx:]
                    else:
                        cols[i] = val[0]
            outf.write('\t'.join(cols) + '\n')
    
    os.replace(temp_file, vcf_file)
    print("Patch applied successfully.")
except Exception as e:
    print(f"Patch failed: {e}")
    if os.path.exists(temp_file):
        os.remove(temp_file)
EOF

# Final index for chrX
bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

echo "Done."
