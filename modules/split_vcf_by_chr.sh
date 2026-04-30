#!/bin/bash
#SBATCH --job-name=split_vcf_by_chr
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

#########################################################################
# This script:
# 1. Extracts the "Pass QC" samples from the PLINK2 filtered dataset.
# 2. Creates a clean male list for those samples.
# 3. Splits the joint VCF by chromosome (1-22, X, Y).
# 4. Forces haploid genotypes for males on ChrX to satisfy TOPMed/eQTL requirements.
#########################################################################

export LC_ALL=C
module load BCFtools/1.21

# Paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"

# 1. Prep lists
awk 'NR>1 {print $1}' $PSAM_FILE > ${OUTPUT_DIR}/pass_ids.txt
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/full_meta.txt
grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="male" {print $1}' > ${OUTPUT_DIR}/males.txt

# 2. Extract ChrX
echo "Extracting chrX..."
bcftools view -r chrX -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# 3. Apply PAR-aware Patch
echo "Applying PAR-aware ploidy patch..."
python3 << 'EOF'
import gzip
import os

vcf_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/target_ALS_chrX.vcf.gz"
male_file = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/males.txt"
temp_file = vcf_file + ".tmp.gz"

# GRCh38 PAR Coordinates (TOPMed / Michigan Server standards)
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
        # Check if position falls within PAR1 or PAR2
        is_par = any(start <= pos <= end for start, end in PAR)
        
        for i in range(9, len(cols)):
            if samples[i-9] in males:
                val = cols[i]
                allele = val[0]
                # Preserve format tags (AD, DP, GQ, etc.)
                sep_idx = val.find(":")
                rest = val[sep_idx:] if sep_idx != -1 else ""
                
                if is_par:
                    # PAR: Force Diploid (e.g., 0/0 or ./.)
                    if "/" not in val[:3] and "|" not in val[:3]:
                        cols[i] = f"{allele}/{allele}{rest}"
                else:
                    # Non-PAR: Force Haploid (e.g., 0 or .)
                    cols[i] = f"{allele}{rest}"
        outf.write('\t'.join(cols) + '\n')

os.replace(temp_file, vcf_file)
EOF

# 4. Re-compress to BGZF and Index (Required for server compatibility)
echo "Fixing compression and indexing..."
bcftools view ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz -Oz -o ${OUTPUT_DIR}/target_ALS_chrX_final.vcf.gz
mv ${OUTPUT_DIR}/target_ALS_chrX_final.vcf.gz ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz
bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# 5. Extract Autosomes
for chr in {1..22} Y; do
    echo "Processing chr${chr}..."
    bcftools view -r chr${chr} -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

echo "Done."
