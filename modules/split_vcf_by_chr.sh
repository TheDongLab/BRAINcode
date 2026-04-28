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

# --- STEP 1: GENERATE SEX MAPPING ---
echo "Generating sex mapping..."
awk -F',' 'NR>1 {print $2, tolower($5)}' $METADATA | sort | uniq | while read id sex; do
    if [[ "$sex" == "male" ]]; then
        echo "$id 1"  # 1 = Male for bcftools ploidy
    elif [[ "$sex" == "female" ]]; then
        echo "$id 2"  # 2 = Female for bcftools ploidy
    fi
done > ${OUTPUT_DIR}/sex_map.txt

# Create a strict ploidy definition file
# Format: Chromosome, Start, End, Sex, Ploidy
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
        
        # Step A: Extract just chrX first
        bcftools view -r chrX $INPUT_VCF -Oz -o ${OUTPUT_DIR}/temp_chrX.vcf.gz
        bcftools index -t ${OUTPUT_DIR}/temp_chrX.vcf.gz
        
        # Step B: Apply fixploidy to the extracted file
        # The '-m f' forces the conversion based on the ploidy file
        bcftools +fixploidy ${OUTPUT_DIR}/temp_chrX.vcf.gz -Oz -o $OUT_VCF -- -p ${OUTPUT_DIR}/ploidy_definition.txt -s ${OUTPUT_DIR}/sex_map.txt
        
        rm ${OUTPUT_DIR}/temp_chrX.vcf.gz ${OUTPUT_DIR}/temp_chrX.vcf.gz.tbi
    else
        bcftools view -r chr${chr} $INPUT_VCF -Oz -o $OUT_VCF
    fi
    
    # Using -f to force overwrite existing index files to avoid the [E::main_vcfindex] error
    bcftools index -f -t $OUT_VCF
done

echo "Process complete. Upload the NEW target_ALS_chrX.vcf.gz."
