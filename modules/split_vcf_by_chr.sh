#!/bin/bash
#SBATCH --job-name=split_vcf_TargetALS
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00

module load BCFtools

# Define paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"

mkdir -p $OUTPUT_DIR

# --- STEP 1: GENERATE SEX MAPPING ---
echo "Generating sex mapping from TargetALS metadata using Subject IDs..."
# Column 2: externalsubjectid (Matches VCF), Column 5: Sex.
awk -F',' 'NR>1 {print $2, tolower($5)}' $METADATA | sort | uniq | while read id sex; do
    if [[ "$sex" == "male" ]]; then
        echo "$id M"
    elif [[ "$sex" == "female" ]]; then
        echo "$id F"
    fi
done > ${OUTPUT_DIR}/sex_map.txt

# Create temporary list of male samples for setGT
grep " M$" ${OUTPUT_DIR}/sex_map.txt | cut -d' ' -f1 > ${OUTPUT_DIR}/males.txt

# --- STEP 2: LOOP THROUGH CHROMOSOMES ---
for chr in {1..22} X Y; do
    echo "Processing chromosome: chr${chr}..."
    OUT_VCF="${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz"
    
    if [ "$chr" == "X" ]; then
        echo "Applying forced haploidization to chrX for male samples..."
        # Extract chrX and immediately force male samples to haploid (strips the / or |)
        bcftools view -r chrX $INPUT_VCF -Ou | \
        bcftools +setGT -- -t q -n h -s ${OUTPUT_DIR}/males.txt | \
        bcftools view -Oz -o $OUT_VCF
    else
        # Standard extraction for autosomes and Y
        bcftools view -r chr${chr} $INPUT_VCF -Oz -o $OUT_VCF
    fi
    
    bcftools index -t $OUT_VCF
done

# Cleanup temporary file
rm ${OUTPUT_DIR}/males.txt

echo "Process complete. Files are in: $OUTPUT_DIR"
echo "You only need to upload the NEW target_ALS_chrX.vcf.gz to the server."
