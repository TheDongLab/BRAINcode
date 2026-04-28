#!/bin/bash
#SBATCH --job-name=split_vcf_TargetALS
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00

module load BCFtools

# Define paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"

mkdir -p $OUTPUT_DIR

# --- STEP 1: GENERATE SEX MAPPING (with Dash -> Underscore fix) ---
echo "Generating sex mapping from TargetALS metadata..."
awk -F',' 'NR>1 {print $1, tolower($5)}' $METADATA | sed 's/-/_/g' | while read id sex; do
    if [[ "$sex" == "male" ]]; then
        echo "$id M"
    elif [[ "$sex" == "female" ]]; then
        echo "$id F"
    fi
done > ${OUTPUT_DIR}/sex_map.txt

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
        echo "Applying ploidy fix to chrX using TargetALS sex mapping..."
        # Pipe directly through fixploidy to resolve ambiguous 0/1 calls in males
        bcftools view -r chrX $INPUT_VCF -Ou | \
        bcftools +fixploidy -- -p ${OUTPUT_DIR}/ploidy_rules.txt -s ${OUTPUT_DIR}/sex_map.txt | \
        bcftools view -Oz -o $OUT_VCF
    else
        # Standard extraction for autosomes and Y
        bcftools view -r chr${chr} $INPUT_VCF -Oz -o $OUT_VCF
    fi
    
    bcftools index -t $OUT_VCF
done

echo "Splitting and ChrX correction complete. Check ${OUTPUT_DIR} for results."
