#!/bin/bash
#SBATCH --job-name=split_vcf_TargetALS
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00

#####################
# Run this script AFTER genotyping QC (i.e, after convert_vcf_to_plink.sh)
#####################

module load BCFtools

# Define paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"

mkdir -p $OUTPUT_DIR

# --- STEP 1: GENERATE SEX MAPPING ---
# tr -d '\r' removes hidden Windows characters
# awk maps Column 2 (ID) and Column 5 (Sex)
echo "Generating sex mapping..."
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | \
sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | sort | uniq | while read id sex; do
    if [[ "$sex" == "male" ]]; then echo "$id 1"; elif [[ "$sex" == "female" ]]; then echo "$id 2"; fi
done > ${OUTPUT_DIR}/sex_map.txt

# Create the ploidy definition file for hg38
cat <<EOF > ${OUTPUT_DIR}/ploidy_definition.txt
X 1 155701382 M 1
X 1 155701382 F 2
EOF

# --- STEP 2: LOOP THROUGH CHROMOSOMES ---
for chr in {1..22} X Y; do
    echo "Processing chromosome: chr${chr}..."
    OUT_VCF="${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz"
    
    if [ "$chr" == "X" ]; then
        echo "Applying ploidy fix to chrX..."
        # Extract and fix ploidy in one pipe to save disk space/I/O
        bcftools view -r chrX $INPUT_VCF -Ou | \
        bcftools +fixploidy -- -p ${OUTPUT_DIR}/ploidy_definition.txt -s ${OUTPUT_DIR}/sex_map.txt | \
        bcftools view -Oz -o $OUT_VCF
    else
        # Standard extraction for autosomes and Y
        bcftools view -r chr${chr} $INPUT_VCF -Oz -o $OUT_VCF
    fi
    
    # -f forces index overwrite to avoid the existing index error
    bcftools index -f -t $OUT_VCF
done

echo "Process complete. Files are in: $OUTPUT_DIR"
