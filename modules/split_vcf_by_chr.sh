#!/bin/bash
#SBATCH --job-name=split_vcf_TargetALS
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00

#####################
# Run this script AFTER genotyping QC
#####################

module load BCFtools

# Define paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"

mkdir -p $OUTPUT_DIR

# --- STEP 1: GENERATE CLEAN MALE LIST ---
# We only care about males because they are the ones causing "ambiguous" errors on ChrX
echo "Generating male sample list..."
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" && (tolower($5) == "male") {print $2}' | \
sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | sort | uniq > ${OUTPUT_DIR}/males.txt

# --- STEP 2: LOOP THROUGH CHROMOSOMES ---
for chr in {1..22} X Y; do
    echo "Processing chromosome: chr${chr}..."
    OUT_VCF="${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz"
    
    if [ "$chr" == "X" ]; then
        echo "Applying definitive haploid fix to chrX..."
        
        # 1. Get header
        bcftools view -r chrX $INPUT_VCF -h > ${OUTPUT_DIR}/header.txt
        
        # 2. Re-code male genotypes to be haploid (strip / and |)
        # This awk script finds which columns are male and removes the slashes from their genotypes
        bcftools view -r chrX $INPUT_VCF -H | \
        awk -v males="$(cat ${OUTPUT_DIR}/males.txt | tr '\n' ',')" -v head="$(grep "^#CHROM" ${OUTPUT_DIR}/header.txt)" '
        BEGIN {
            split(males, m_arr, ",");
            for (i in m_arr) is_male[m_arr[i]] = 1;
            split(head, h_cols, "\t");
            for (i=10; i<=length(h_cols); i++) if (h_cols[i] in is_male) male_col[i] = 1;
            FS=OFS="\t";
        }
        {
            for (i=10; i<=NF; i++) {
                if (male_col[i]) gsub(/[\/|]/, "", $i);
            }
            print $0;
        }' | cat ${OUTPUT_DIR}/header.txt - | bcftools view -Oz -o $OUT_VCF
        
        rm ${OUTPUT_DIR}/header.txt
    else
        # Standard extraction for autosomes and Y
        bcftools view -r chr${chr} $INPUT_VCF -Oz -o $OUT_VCF
    fi
    
    # -f forces index overwrite
    bcftools index -f -t $OUT_VCF
done

echo "Process complete. Files are in: $OUTPUT_DIR"
