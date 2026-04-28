#!/bin/bash
#SBATCH --job-name=split_vcf_QC_filtered
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

module load BCFtools

# Define paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"

mkdir -p $OUTPUT_DIR

# --- STEP 1: EXTRACT SAMPLES ---
PASS_SAMPLES="${OUTPUT_DIR}/samples_pass_QC.txt"
awk 'NR>1 {print $1}' $PSAM_FILE > $PASS_SAMPLES

# --- STEP 2: GENERATE MALE LIST ---
MALE_LIST="${OUTPUT_DIR}/males_pass_qc.txt"
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" && (tolower($5) == "male") {print $2}' | \
sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | grep -Fwf $PASS_SAMPLES | sort | uniq > $MALE_LIST

# --- STEP 3: LOOP THROUGH CHROMOSOMES ---
for chr in {1..22} X Y; do
    echo "Processing chromosome: chr${chr}..."
    OUT_VCF="${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz"
    
    if [ "$chr" == "X" ]; then
        echo "Applying first-allele haploid fix to chrX..."
        bcftools view -r chrX -S $PASS_SAMPLES --force-samples $INPUT_VCF -h > ${OUTPUT_DIR}/temp_header.txt
        
        bcftools view -r chrX -S $PASS_SAMPLES --force-samples $INPUT_VCF -H | \
        awk -v males="$(cat $MALE_LIST | tr '\n' ',')" -v head="$(grep "^#CHROM" ${OUTPUT_DIR}/temp_header.txt)" '
        BEGIN {
            split(males, m_arr, ",");
            for (i in m_arr) is_male[m_arr[i]] = 1;
            split(head, h_cols, "\t");
            for (i=10; i<=length(h_cols); i++) if (h_cols[i] in is_male) male_col[i] = 1;
            FS=OFS="\t";
        }
        {
            for (i=10; i<=NF; i++) {
                if (male_col[i]) {
                    # Take only the very first character (0, 1, or .)
                    # This ensures 0/1 -> 0, 1/1 -> 1, and ./. -> .
                    $i = substr($i, 1, 1);
                }
            }
            print $0;
        }' | cat ${OUTPUT_DIR}/temp_header.txt - | bcftools view -Oz -o $OUT_VCF
        
        rm ${OUTPUT_DIR}/temp_header.txt
    else
        bcftools view -r chr${chr} -S $PASS_SAMPLES --force-samples $INPUT_VCF -Oz -o $OUT_VCF
    fi
    
    bcftools index -f -t $OUT_VCF
done
