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
module load HTSlib/1.21
module load PLINK2

############################################
# INPUTS
############################################

VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
REF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"
USER_DATA_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
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
-o "${PREFIX}.step1.bcf"

bcftools index -f "${PREFIX}.step1.bcf"

############################################
# STEP 4: SPLIT MULTIALLELICS + LEFT NORMALIZE
############################################

echo "==========================================="
echo "STEP 4: SPLITTING + NORMALIZING"
echo "==========================================="

bcftools norm \
-f "${REF}" \
-m -any \
"${PREFIX}.step1.bcf" \
-Ou \
-o "${PREFIX}.step2.bcf"

bcftools index -f "${PREFIX}.step2.bcf"

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
-o "${PREFIX}.step3.bcf"

bcftools index -f "${PREFIX}.step3.bcf"

############################################
# STEP 6: REMOVE LARGE/COMPLEX INDELS
############################################
#
# PLINK2 can fail on huge alleles and highly
# repetitive complex indels.
#
# Keep:
#   SNPs
#   small/medium indels
#
# Remove:
#   alleles > 50bp
#
############################################

echo "==========================================="
echo "STEP 6: REMOVING LARGE/COMPLEX INDELS"
echo "==========================================="

bcftools filter \
-i '
(
    strlen(REF)=1 &&
    strlen(ALT)=1
)
||
(
    length(REF)<=50 &&
    length(ALT)<=50
)
' \
"${PREFIX}.step3.bcf" \
-Ob \
-o "${PREFIX}.step4.bcf"

bcftools index -f "${PREFIX}.step4.bcf"

############################################
# STEP 7: REMOVE MONOMORPHIC SITES
############################################

echo "==========================================="
echo "STEP 7: REMOVING MONOMORPHIC SITES"
echo "==========================================="

bcftools view \
-c 1:minor \
"${PREFIX}.step4.bcf" \
-Ou \
-o "${PREFIX}.step5.bcf"

bcftools index -f "${PREFIX}.step5.bcf"

############################################
# STEP 8: CONVERT '.' TO './.'
############################################

echo "==========================================="
echo "STEP 8: FIXING MISSING GT"
echo "==========================================="

bcftools +setGT \
"${PREFIX}.step5.bcf" \
-Ou \
-- \
-t q \
-n "./." \
-i 'GT="."' \
> "${PREFIX}.step6.bcf"

bcftools index -f "${PREFIX}.step6.bcf"

############################################
# STEP 9: REMOVE PHASING TAGS
############################################

echo "==========================================="
echo "STEP 9: REMOVING PHASING TAGS"
echo "==========================================="

bcftools annotate \
-x FORMAT/PGT,FORMAT/PID \
"${PREFIX}.step6.bcf" \
-Ou \
-o "${PREFIX}.step7.bcf"

bcftools index -f "${PREFIX}.step7.bcf"

############################################
# STEP 10: WRITE FINAL VCF
############################################

echo "==========================================="
echo "STEP 10: WRITING FINAL VCF"
echo "==========================================="

bcftools view \
"${PREFIX}.step7.bcf" \
-Oz \
-o "${PREFIX}.cleaned.vcf.gz"

tabix -f -p vcf "${PREFIX}.cleaned.vcf.gz"

############################################
# STEP 11: FINAL GT SUMMARY
############################################

echo "==========================================="
echo "STEP 11: FINAL GT SUMMARY"
echo "==========================================="

bcftools query \
-f '[%GT\n]' \
"${PREFIX}.cleaned.vcf.gz" \
| sort \
| uniq -c \
> "${PREFIX}.final_GT_summary.txt"

############################################
# STEP 12: BAD GT TOKEN CHECK
############################################

echo "==========================================="
echo "STEP 12: BAD GT CHECK"
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
# STEP 13: COLUMN CONSISTENCY CHECK
############################################

echo "==========================================="
echo "STEP 13: COLUMN CONSISTENCY CHECK"
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
# STEP 14: SYMBOLIC ALLELE CHECK
############################################

echo "==========================================="
echo "STEP 14: SYMBOLIC ALLELE CHECK"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\n' \
"${PREFIX}.cleaned.vcf.gz" \
| awk '
$4 ~ /</ ||
$4 == "*" ||
$4 ~ /,/ {
    print
}' \
> "${PREFIX}.remaining_symbolic_alleles.txt"

############################################
# STEP 15: LARGE ALLELE CHECK
############################################

echo "==========================================="
echo "STEP 15: LARGE ALLELE CHECK"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\n' \
"${PREFIX}.cleaned.vcf.gz" \
| awk '
length($3)>50 ||
length($4)>50
' \
> "${PREFIX}.remaining_large_alleles.txt"

############################################
# STEP 16: MULTIALLELIC CHECK
############################################

echo "==========================================="
echo "STEP 16: MULTIALLELIC CHECK"
echo "==========================================="

bcftools view \
-H \
"${PREFIX}.cleaned.vcf.gz" \
| awk '
$5 ~ /,/
' \
> "${PREFIX}.remaining_multiallelic.txt"

############################################
# STEP 17: VCF VALIDATION
############################################

echo "==========================================="
echo "STEP 17: VCF VALIDATION"
echo "==========================================="

bcftools view \
"${PREFIX}.cleaned.vcf.gz" \
> /dev/null

############################################
# STEP 18: PLINK2 IMPORT TEST
############################################
#
# IMPORTANT:
# This is PLINK2, not PLINK1.9
#
# We use:
#   --make-pgen
#
############################################

echo "==========================================="
echo "STEP 18: PLINK2 IMPORT TEST"
echo "==========================================="

plink2 \
--vcf "${PREFIX}.cleaned.vcf.gz" \
--make-pgen \
--new-id-max-allele-len 100 truncate \
--snps-only just-acgt \
--split-par b38 \
--out "${PREFIX}.plink2_test"

############################################
# STEP 19: FINAL VARIANT COUNTS
############################################

echo "==========================================="
echo "STEP 19: FINAL COUNTS"
echo "==========================================="

{
echo "RAW_VARIANTS"
bcftools view -H "${PREFIX}.female.raw.vcf.gz" | wc -l

echo "FINAL_VARIANTS"
bcftools view -H "${PREFIX}.cleaned.vcf.gz" | wc -l
} \
> "${PREFIX}.variant_counts.txt"

############################################
# DONE
############################################

echo "==========================================="
echo "DONE"
echo "==========================================="

echo ""
echo "FINAL FILE:"
echo "${PREFIX}.cleaned.vcf.gz"

echo ""
echo "PLINK2 FILES:"
echo "${PREFIX}.plink2_test.pgen"
echo "${PREFIX}.plink2_test.pvar"
echo "${PREFIX}.plink2_test.psam"

echo ""
echo "VALIDATION FILES:"
echo "${PREFIX}.initial_GT_summary.txt"
echo "${PREFIX}.final_GT_summary.txt"
echo "${PREFIX}.bad_GT_tokens.txt"
echo "${PREFIX}.bad_lines.txt"
echo "${PREFIX}.remaining_symbolic_alleles.txt"
echo "${PREFIX}.remaining_large_alleles.txt"
echo "${PREFIX}.remaining_multiallelic.txt"
echo "${PREFIX}.variant_counts.txt"
