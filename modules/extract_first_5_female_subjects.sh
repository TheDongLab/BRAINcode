#!/bin/bash

# --- PATHS ---
VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
FEMALES_SRC="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/females.txt"
REF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"
OUT_VCF="${OUT_DIR}/TargetALS_chrX_first5_Females.vcf.gz"

module load BCFtools/1.21

# 1. Create a temporary list of the first 5 female IDs
head -n 5 "$FEMALES_SRC" > "${OUT_DIR}/top5_female_ids.txt"

# 2. Extract, Normalize, and Zip
# We include 'chrX' to limit the scope and 'norm' to ensure the VCF is clean
bcftools view -S "${OUT_DIR}/top5_female_ids.txt" "$VCF_IN" chrX | \
bcftools norm -m -any -f "$REF" | \
bcftools sort -Oz -o "$OUT_VCF"

# 3. Index the resulting VCF
bcftools index -f -t "$OUT_VCF"

# Verification
echo "VCF Created: $OUT_VCF"
echo "Samples included:"
bcftools query -l "$OUT_VCF"
