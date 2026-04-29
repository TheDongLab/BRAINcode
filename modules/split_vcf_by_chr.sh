#!/bin/bash
#SBATCH --job-name=split_vcf_by_chr
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

#########################################################################
# This script:
# 1. Extracts the "Pass QC" samples from the PLINK2 filtered dataset.
# 2. Creates a clean male list for those samples.
# 3. Splits the joint VCF by chromosome (1-22, X, Y).
# 4. Forces haploid genotypes for males on ChrX to satisfy TOPMed/eQTL requirements.
#########################################################################

module load BCFtools/1.21

# Paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"

mkdir -p $OUTPUT_DIR

# --- STEP 1: PREP REQUISITES ---
echo "Creating sex-specific sample lists..."

# Extract IDs from PLINK pass list
RAW_PASS="${OUTPUT_DIR}/raw_pass_ids.txt"
awk 'NR>1 {print $1}' $PSAM_FILE > $RAW_PASS

# Create Female and Male lists
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/full_meta.txt
grep -Fwf $RAW_PASS ${OUTPUT_DIR}/full_meta.txt | awk '$2=="female" {print $1}' > ${OUTPUT_DIR}/females.txt
grep -Fwf $RAW_PASS ${OUTPUT_DIR}/full_meta.txt | awk '$2=="male" {print $1}' > ${OUTPUT_DIR}/males.txt

# Final sample list for subsetting
PASS_SAMPLES="${OUTPUT_DIR}/pass_samples.txt"
cat ${OUTPUT_DIR}/females.txt ${OUTPUT_DIR}/males.txt > $PASS_SAMPLES

# --- STEP 2: CHRX (PYTHON RE-WRITE) ---
echo "Processing chrX: Subsetting and then applying Python brute-force fix..."

# 2a. Initial subset (likely still has inconsistent ploidy)
bcftools view -r chrX -S $PASS_SAMPLES --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/tmp_chrX_prefixed.vcf.gz

# 2b. Python brute-force ploidy correction
# This reads the VCF and FORCES slashes for females and NO slashes for males.
python3 << 'EOF'
import gzip
import os

output_dir = "/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
female_file = os.path.join(output_dir, "females.txt")
vcf_in = os.path.join(output_dir, "tmp_chrX_prefixed.vcf.gz")
vcf_out = os.path.join(output_dir, "target_ALS_chrX.vcf.gz")

with open(female_file, 'r') as f:
    females = set(line.strip() for line in f)

with gzip.open(vcf_in, 'rt') as inf, gzip.open(vcf_out, 'wt') as outf:
    for line in inf:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
            outf.write(line)
            continue
        
        cols = line.strip().split('\t')
        info_fields = cols[:9]
        genotypes = cols[9:]
        new_genotypes = []
        
        for i, gt_field in enumerate(genotypes):
            sample_name = samples[i]
            parts = gt_field.split(':')
            # Extract only the alleles, remove any current separators
            raw_alleles = parts[0].replace('/', '').replace('|', '')
            
            if sample_name in females:
                # Force Diploid
                if raw_alleles == ".":
                    parts[0] = "./."
                elif len(raw_alleles) == 1:
                    parts[0] = f"{raw_alleles}/{raw_alleles}"
                else:
                    parts[0] = f"{raw_alleles[0]}/{raw_alleles[1]}"
            else:
                # Force Haploid
                parts[0] = raw_alleles[0] if raw_alleles else "."
            
            new_genotypes.append(':'.join(parts))
        
        outf.write('\t'.join(info_fields) + '\t' + '\t'.join(new_genotypes) + '\n')
EOF

bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# --- STEP 3: ALL OTHER CHROMOSOMES ---
for chr in {1..22} Y; do
    echo "Processing chr${chr}..."
    bcftools view -r chr${chr} -S $PASS_SAMPLES --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

# --- STEP 4: CLEANUP ---
rm -f ${OUTPUT_DIR}/tmp_chrX_prefixed.vcf.gz ${OUTPUT_DIR}/raw_pass_ids.txt ${OUTPUT_DIR}/full_meta.txt
echo "Success. Pipeline complete."
