#!/bin/bash
#SBATCH --job-name=clean_chrX_topmed
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/clean_chrX_topmed.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/clean_chrX_topmed.err

set -euo pipefail

############################################
# MODULES
############################################

module purge
module load BCFtools/1.21
module load PLINK/1.9b_7.11-x86_64

############################################
# INPUTS
############################################

VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
PSAM_IN="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"
REF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"
USER_DATA_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
MALES_SRC="${USER_DATA_DIR}/males.txt"
FEMALES_SRC="${USER_DATA_DIR}/females.txt"
PREFIX="${OUT_DIR}/target_ALS_chrX_TOPMed"

mkdir -p "${OUT_DIR}"

############################################
# STEP 1: EXTRACT FEMALE chrX
############################################

echo "==========================================="
echo "STEP 1: EXTRACTING FEMALE chrX"
echo "==========================================="

bcftools view \
-S "${FEMALES_SRC}" \
-r chrX,X \
"${VCF_IN}" \
-Oz \
-o "${PREFIX}.female.raw.vcf.gz"

tabix -f -p vcf "${PREFIX}.female.raw.vcf.gz"

############################################
# STEP 2: INITIAL GT SUMMARY
############################################

echo "==========================================="
echo "STEP 2: INITIAL GT SUMMARY"
echo "==========================================="

bcftools query \
-f '[%GT\n]' \
"${PREFIX}.female.raw.vcf.gz" \
| sort \
| uniq -c \
> "${PREFIX}.initial_GT_summary.txt"

############################################
# STEP 3: REMOVE SYMBOLIC / NON_REF ALLELES
############################################

echo "==========================================="
echo "STEP 3: REMOVING SYMBOLIC/NON_REF"
echo "==========================================="

bcftools view \
-e 'ALT="<NON_REF>" || ALT="*" || ALT~"<"' \
"${PREFIX}.female.raw.vcf.gz" \
-Ou \
> "${PREFIX}.step1.bcf"

############################################
# STEP 4: SPLIT MULTIALLELICS
############################################

echo "==========================================="
echo "STEP 4: SPLITTING MULTIALLELICS"
echo "==========================================="

bcftools norm \
-f "${REF}" \
-m -any \
"${PREFIX}.step1.bcf" \
-Ou \
> "${PREFIX}.step2.bcf"

############################################
# STEP 5: KEEP ONLY BIALLELIC SNP/INDEL
############################################

echo "==========================================="
echo "STEP 5: KEEPING BIALLELIC SNP/INDEL"
echo "==========================================="

bcftools view \
-m2 \
-M2 \
-v snps,indels \
"${PREFIX}.step2.bcf" \
-Ou \
> "${PREFIX}.step3.bcf"

############################################
# STEP 6: REMOVE MONOMORPHIC SITES
############################################

echo "==========================================="
echo "STEP 6: REMOVING MONOMORPHIC SITES"
echo "==========================================="

bcftools view \
-c 1:minor \
"${PREFIX}.step3.bcf" \
-Ou \
> "${PREFIX}.step4.bcf"

############################################
# STEP 7: CONVERT '.' TO './.'
############################################

echo "==========================================="
echo "STEP 7: FIXING MISSING GT"
echo "==========================================="

bcftools +setGT \
"${PREFIX}.step4.bcf" \
-Ou \
-- \
-n "./." \
-i 'GT="."' \
> "${PREFIX}.step5.bcf"

############################################
# STEP 8: REMOVE PHASING
############################################

echo "==========================================="
echo "STEP 8: REMOVING PHASING"
echo "==========================================="

bcftools query \
-l \
"${PREFIX}.step5.bcf" \
> "${PREFIX}.samples.txt"

bcftools reheader \
-s "${PREFIX}.samples.txt" \
"${PREFIX}.step5.bcf" \
-Ou \
| bcftools annotate \
-x FORMAT/PGT,FORMAT/PID \
-Ou \
| bcftools +setGT \
-Ou \
-- \
-t q \
-n u \
> "${PREFIX}.step6.bcf"

############################################
# STEP 9: WRITE FINAL VCF
############################################

echo "==========================================="
echo "STEP 9: WRITING FINAL VCF"
echo "==========================================="

bcftools view \
"${PREFIX}.step6.bcf" \
-Oz \
-o "${PREFIX}.cleaned.vcf.gz"

tabix -f -p vcf "${PREFIX}.cleaned.vcf.gz"

############################################
# STEP 10: FINAL GT SUMMARY
############################################

echo "==========================================="
echo "STEP 10: FINAL GT SUMMARY"
echo "==========================================="

bcftools query \
-f '[%GT\n]' \
"${PREFIX}.cleaned.vcf.gz" \
| sort \
| uniq -c \
> "${PREFIX}.final_GT_summary.txt"

############################################
# STEP 11: BAD GT CHECK
############################################

echo "==========================================="
echo "STEP 11: BAD GT CHECK"
echo "==========================================="

bcftools query \
-f '[%SAMPLE\t%GT\n]' \
"${PREFIX}.cleaned.vcf.gz" \
| awk '
$2 != "." &&
$2 != "./." &&
$2 !~ /^[0-9]+\/[0-9]+$/ {
    print
}' \
> "${PREFIX}.bad_GT_tokens.txt"

############################################
# STEP 12: COLUMN CONSISTENCY CHECK
############################################

echo "==========================================="
echo "STEP 12: COLUMN CONSISTENCY CHECK"
echo "==========================================="

zcat "${PREFIX}.cleaned.vcf.gz" \
| awk '
BEGIN{FS="\t"}
/^#CHROM/{
    expected=NF
    next
}
!/^#/ && NF!=expected{
    print "BAD_LINE", NR, NF, expected, $1, $2
}' \
> "${PREFIX}.bad_lines.txt"

############################################
# STEP 13: CHECK FOR REMAINING SYMBOLIC ALLELES
############################################

echo "==========================================="
echo "STEP 13: SYMBOLIC ALLELE CHECK"
echo "==========================================="

bcftools view \
-H \
"${PREFIX}.cleaned.vcf.gz" \
| awk '
$5 ~ /</ || $5 == "*" {
    print
}' \
> "${PREFIX}.remaining_symbolic_alleles.txt"

############################################
# STEP 14: FINAL SUMMARY
############################################

echo "==========================================="
echo "DONE"
echo "==========================================="

echo "FINAL FILE:"
echo "${PREFIX}.cleaned.vcf.gz"

echo ""
echo "VALIDATION FILES:"
echo "${PREFIX}.initial_GT_summary.txt"
echo "${PREFIX}.final_GT_summary.txt"
echo "${PREFIX}.bad_GT_tokens.txt"
echo "${PREFIX}.bad_lines.txt"
echo "${PREFIX}.remaining_symbolic_alleles.txt"
