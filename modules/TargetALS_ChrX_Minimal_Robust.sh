#!/bin/bash
#SBATCH --job-name=TargetALS_ChrX_Minimal_Robust
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/chrX_robust_fix.out

# --- PATHS ---
VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
REF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"
USER_DATA_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"

module load BCFtools/1.21

for SEX in males females; do
    echo "Processing ${SEX^^}..."
    SRC_LIST="${USER_DATA_DIR}/${SEX}.txt"
    PREFIX="${OUT_DIR}/target_ALS_chrX_${SEX^^}"

    # 1. Subset, Remove Symbolics, and Normalize
    # This remains the most crucial "Structural" fix for imputation.
    bcftools view -S "$SRC_LIST" -r chrX "$VCF_IN" \
        -e 'ALT="<NON_REF>" || ALT="*" || ALT~"<"' -Ou | \
    bcftools norm -m -any -f "$REF" -Ou | \
    \
    # 2. Robust Genotype Cleaning
    # Instead of +setGT, we use 'set-GTs' which is built into 'view'.
    # For Males: We zero out HET calls in non-PAR regions.
    # For Females: We pass through.
    if [ "$SEX" == "males" ]; then
        # hg38 non-PAR filter: set genotypes to missing (.) for HETs outside PAR1/PAR2
        bcftools view -e 'GT="het" && (POS<10001 || (POS>2781479 && POS<155701383) || POS>156030895)' -Ou
    else
        cat
    fi | \
    \
    # 3. Format and Zip
    # Removing phasing blocks and forced unphasing (sed)
    bcftools annotate -x FORMAT/PID,FORMAT/PGT -Ou | \
    sed 's/|/\//g' | \
    bcftools view -Oz -o "${PREFIX}.cleaned.vcf.gz"

    tabix -f -p vcf "${PREFIX}.cleaned.vcf.gz"
done
