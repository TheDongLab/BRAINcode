#!/bin/bash
#SBATCH --job-name=split_vcf_fixX
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00

module load BCFtools

# Define paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
# Path to your existing male-stratified expression or metadata for ID extraction
MALE_METADATA="/home/zw529/donglab/data/answer_ALS/consolidated_metadata_full.tsv"

mkdir -p $OUTPUT_DIR

# --- STEP 1: GENERATE PLOIDY FILES FOR CHRX ---
# Create a sample-to-sex mapping from your metadata
# Assuming Column 1 is SampleID and Column 12 is Sex (adjust 'cut' as needed)
echo "Generating sex mapping from metadata..."
grep -i "Male" $MALE_METADATA | cut -f1 | awk '{print $1" M"}' > ${OUTPUT_DIR}/sex_map.txt
grep -i "Female" $MALE_METADATA | cut -f1 | awk '{print $1" F"}' >> ${OUTPUT_DIR}/sex_map.txt

# Create the ploidy rules (hg38 coordinates for X)
cat <<EOF > ${OUTPUT_DIR}/ploidy_rules.txt
X 1 155701382 M 1
X 1 155701382 F 2
EOF

# --- STEP 2: LOOP THROUGH CHROMOSOMES ---
for chr in {1..22} X Y; do
    echo "Processing chromosome: chr${chr}..."
    OUT_VCF="${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz"
    
    if [ "$chr" == "X" ]; then
        echo "Applying ploidy fix to chrX..."
        # Extract, then pipe into fixploidy plugin to resolve 0/1 calls in males
        bcftools view -r chrX $INPUT_VCF -Ou | \
        bcftools +fixploidy -- -p ${OUTPUT_DIR}/ploidy_rules.txt -s ${OUTPUT_DIR}/sex_map.txt | \
        bcftools view -Oz -o $OUT_VCF
    else
        # Standard extraction for autosomes and Y
        bcftools view -r chr${chr} $INPUT_VCF -Oz -o $OUT_VCF
    fi
    
    bcftools index -t $OUT_VCF
done

echo "Splitting and ChrX correction complete. Files in: $OUTPUT_DIR"
