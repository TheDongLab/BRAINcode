#!/bin/bash
#SBATCH --job-name=TargetALS_ChrX_Stratified_TOPMed
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/clean_chrX_stratified_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/clean_chrX_stratified_%j.err

set -euo pipefail

############################################
# MODULES
############################################
module purge
module load BCFtools/1.21
module load HTSlib/1.21
module load PLINK

############################################
# INPUTS
############################################
VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
REF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"
USER_DATA_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"

mkdir -p "${OUT_DIR}"

# START LOOP
for SEX in males females; do
    echo "==========================================="
    echo "STARTING PIPELINE FOR: ${SEX^^}"
    echo "==========================================="

    SRC_LIST="${USER_DATA_DIR}/${SEX}.txt"
    PREFIX="${OUT_DIR}/target_ALS_chrX_${SEX^^}_TOPMed"

    ############################################
    # STEP 1: EXTRACT chrX
    ############################################
    bcftools view -S "${SRC_LIST}" -r chrX,X "${VCF_IN}" -Oz -o "${PREFIX}.raw.vcf.gz"
    tabix -f -p vcf "${PREFIX}.raw.vcf.gz"

    ############################################
    # STEPS 2-4: INITIAL DIAGNOSTICS
    ############################################
    bcftools view -H "${PREFIX}.raw.vcf.gz" | wc -l > "${PREFIX}.initial_variant_count.txt"
    bcftools stats "${PREFIX}.raw.vcf.gz" > "${PREFIX}.initial.stats.txt"

    ############################################
    # STEP 5: REMOVE SYMBOLIC / NON_REF ALLELES
    ############################################
    bcftools view -e 'ALT="<NON_REF>" || ALT="*" || ALT~"<"' "${PREFIX}.raw.vcf.gz" -Ou > "${PREFIX}.step1.bcf"

    ############################################
    # STEP 6: SPLIT MULTIALLELICS + LEFT NORMALIZE
    ############################################
    bcftools norm -f "${REF}" -m -any "${PREFIX}.step1.bcf" -Ou > "${PREFIX}.step2.bcf"

    ############################################
    # STEP 8: KEEP ONLY BIALLELIC SNP/INDEL
    ############################################
    bcftools view -m2 -M2 -v snps,indels "${PREFIX}.step2.bcf" -Ou > "${PREFIX}.step3.bcf"

    ############################################
    # STEP 9: REMOVE VERY LARGE INDELS
    ############################################
    bcftools view -i 'strlen(REF)<=50 && strlen(ALT)<=50' "${PREFIX}.step3.bcf" -Ou > "${PREFIX}.step4.bcf"

    ############################################
    # STEP 10: REMOVE MONOMORPHIC SITES
    ############################################
    bcftools view -c 1:minor "${PREFIX}.step4.bcf" -Ou > "${PREFIX}.step5.bcf"

    ############################################
    # STEP 10.5: MALE HEMIZYGOUS CLEANING (NEW)
    ############################################
    if [ "$SEX" == "males" ]; then
        echo "Applying Male-specific non-PAR HET filter..."
        bcftools view -e 'GT="het" && (POS<10001 || (POS>2781479 && POS<155701383) || POS>156030895)' \
        "${PREFIX}.step5.bcf" -Ou > "${PREFIX}.step5_5.bcf"
    else
        echo "Passing Females through Step 10.5 without modification..."
        cp "${PREFIX}.step5.bcf" "${PREFIX}.step5_5.bcf"
    fi

    ############################################
    # STEP 11: FIX MISSING GT TOKENS
    ############################################
    bcftools +setGT "${PREFIX}.step5_5.bcf" -Ou -- -t q -n "./." -i 'GT="."' > "${PREFIX}.step6.bcf"

    ############################################
    # STEP 12: REMOVE PHASING TAGS
    ############################################
    bcftools annotate -x FORMAT/PID,FORMAT/PGT "${PREFIX}.step6.bcf" -Ou > "${PREFIX}.step7.bcf"

    ############################################
    # STEP 13: CONVERT PHASED GT TO UNPHASED
    ############################################
    bcftools view "${PREFIX}.step7.bcf" | sed 's/|/\//g' | bcftools view -Ou > "${PREFIX}.step8.bcf"

    ############################################
    # STEP 14: FINAL VCF
    ############################################
    bcftools view "${PREFIX}.step8.bcf" -Oz -o "${PREFIX}.cleaned.vcf.gz"
    tabix -f -p vcf "${PREFIX}.cleaned.vcf.gz"

    ############################################
    # STEPS 15-21: FINAL QUALITY CHECKS
    ############################################
    bcftools query -f '[%GT\n]' "${PREFIX}.cleaned.vcf.gz" | sort | uniq -c > "${PREFIX}.final_GT_summary.txt"
    bcftools view -H "${PREFIX}.cleaned.vcf.gz" | wc -l > "${PREFIX}.final_variant_count.txt"
    
    # Check for remaining bad tokens or symbolic alleles
    bcftools query -f '[%SAMPLE\t%GT\n]' "${PREFIX}.cleaned.vcf.gz" | awk '$2 != "." && $2 != "./." && $2 !~ /^[0-9]+\/[0-9]+$/ {print}' > "${PREFIX}.bad_GT_tokens.txt"
    bcftools view -H "${PREFIX}.cleaned.vcf.gz" | awk '$5 ~ /</ || $5 == "*" {print}' > "${PREFIX}.remaining_symbolic_alleles.txt"

    ############################################
    # STEP 22: PLINK VALIDATION
    ############################################
    plink --vcf "${PREFIX}.cleaned.vcf.gz" --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --out "${PREFIX}.plink_test"

    echo "COMPLETED: ${SEX^^}"
done

echo "==========================================="
echo "ALL COHORTS PROCESSED SUCCESSFULLY"
echo "==========================================="
