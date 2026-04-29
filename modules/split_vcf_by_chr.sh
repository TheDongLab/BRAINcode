#!/bin/bash
#SBATCH --job-name=split_vcf_by_chr
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

# Exit if any command fails
set -e

# Prevent multibyte/locale errors by forcing raw byte handling
export LC_ALL=C

module load BCFtools/1.21

# Paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"

mkdir -p $OUTPUT_DIR

# 1. Prep Sex Lists
echo "Generating sex-specific sample lists..."
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/full_meta.txt
awk 'NR>1 {print $1}' $PSAM_FILE > ${OUTPUT_DIR}/pass_ids.txt

FEMALES=$(grep -Fwf ${OUTPUT_DIR}/pass_ids.txt ${OUTPUT_DIR}/full_meta.txt | awk '$2=="female" {print $1}' | xargs)

# 2. Process ChrX
echo "Applying brute-force string replacement on chrX..."
# Use -Ou for uncompressed output to speed up the pipe to awk
bcftools view -r chrX -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Ou | \
awk -v fem="$FEMALES" '
BEGIN {
    FS="\t"; OFS="\t";
    n=split(fem, f_arr, " "); for(i=1; i<=n; i++) is_fem[f_arr[i]]=1;
}
/^##/ { print; next }
/^#CHROM/ {
    print;
    for(i=9; i<=NF; i++) sample_map[i]=$i;
    next
}
{
    for(i=10; i<=NF; i++) {
        # Grab the first character (the allele)
        allele = substr($i, 1, 1);
        
        # Capture any format tags that follow (AD, DP, GQ, etc.)
        idx = index($i, ":");
        rest = (idx > 0) ? substr($i, idx) : "";

        if (sample_map[i] in is_fem) {
            # Force Diploid string: ./. or 0/0 or 1/1
            new_gt = (allele == ".") ? "./." : allele "/" allele;
            $i = new_gt rest;
        } else {
            # Force Haploid string: . or 0 or 1
            $i = allele rest;
        }
    }
    print;
}' | bcftools view -Oz -o ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# 3. Autosomes & Y
echo "Processing remaining chromosomes..."
for chr in {1..22} Y; do
    echo "Processing chr${chr}..."
    bcftools view -r chr${chr} -S ${OUTPUT_DIR}/pass_ids.txt --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

echo "SUCCESS. All files processed, fixed, and indexed."
