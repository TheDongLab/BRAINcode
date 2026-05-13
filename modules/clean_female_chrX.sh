#!/bin/bash
#SBATCH --job-name=clean_chrX_topmed
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=clean_chrX_topmed.out
#SBATCH --error=clean_chrX_topmed.err

set -euo pipefail

############################################
# INPUTS
############################################

VCF="target_ALS_chrX.vcf.gz"
REF="reference.fa"
PREFIX="target_ALS_chrX.TOPMed"

############################################
# LOAD MODULES
############################################

module load BCFtools
module load htslib

############################################
# STEP 1: INITIAL SUMMARY
############################################

echo "==========================================="
echo "INITIAL GT SUMMARY"
echo "==========================================="

bcftools query -f '[%GT\n]' "$VCF" \
| sort \
| uniq -c \
> "${PREFIX}.initial_GT_summary.txt"

############################################
# STEP 2: REMOVE SYMBOLIC / GVCF ALLELES
############################################

echo "==========================================="
echo "REMOVING NON_REF AND SYMBOLIC ALLELES"
echo "==========================================="

bcftools view \
-e 'ALT="<NON_REF>" || ALT="*" || ALT~"<"' \
"$VCF" \
-Ou \
> "${PREFIX}.step1.bcf"

############################################
# STEP 3: SPLIT MULTIALLELICS
############################################

echo "==========================================="
echo "NORMALIZING MULTIALLELICS"
echo "==========================================="

bcftools norm \
-f "$REF" \
-m -any \
"${PREFIX}.step1.bcf" \
-Ou \
> "${PREFIX}.step2.bcf"

############################################
# STEP 4: KEEP ONLY BIALLELIC VARIANTS
############################################

echo "==========================================="
echo "KEEPING ONLY BIALLELIC SNPs/INDELs"
echo "==========================================="

bcftools view \
-m2 \
-M2 \
-v snps,indels \
"${PREFIX}.step2.bcf" \
-Ou \
> "${PREFIX}.step3.bcf"

############################################
# STEP 5: REMOVE MONOMORPHIC SITES
############################################

echo "==========================================="
echo "REMOVING MONOMORPHIC SITES"
echo "==========================================="

bcftools view \
-c 1:minor \
"${PREFIX}.step3.bcf" \
-Ou \
> "${PREFIX}.step4.bcf"

############################################
# STEP 6: FORCE DIPLOID MISSINGNESS
############################################

echo "==========================================="
echo "FIXING MISSING GT REPRESENTATION"
echo "==========================================="

bcftools +setGT \
"${PREFIX}.step4.bcf" \
-- \
-n "./." \
-i 'GT="."' \
-Ou \
> "${PREFIX}.step5.bcf"

############################################
# STEP 7: REMOVE PHASING
############################################

echo "==========================================="
echo "CONVERTING PHASED TO UNPHASED"
echo "==========================================="

bcftools query -l "${PREFIX}.step5.bcf" > samples.txt

bcftools reheader \
-s samples.txt \
"${PREFIX}.step5.bcf" \
-Ou \
| bcftools annotate \
-x FORMAT/PGT,FORMAT/PID \
-Ou \
| bcftools +setGT \
-- \
-t q \
-n u \
-Ou \
> "${PREFIX}.step6.bcf"

############################################
# STEP 8: WRITE FINAL VCF
############################################

echo "==========================================="
echo "WRITING FINAL VCF"
echo "==========================================="

bcftools view \
"${PREFIX}.step6.bcf" \
-Oz \
-o "${PREFIX}.cleaned.vcf.gz"

tabix -f -p vcf "${PREFIX}.cleaned.vcf.gz"

############################################
# STEP 9: FINAL VALIDATION
############################################

echo "==========================================="
echo "FINAL GT SUMMARY"
echo "==========================================="

bcftools query -f '[%GT\n]' "${PREFIX}.cleaned.vcf.gz" \
| sort \
| uniq -c \
> "${PREFIX}.final_GT_summary.txt"

echo "==========================================="
echo "CHECKING FOR BAD GT TOKENS"
echo "==========================================="

bcftools query -f '[%SAMPLE\t%GT\n]' "${PREFIX}.cleaned.vcf.gz" \
| awk '
$2 != "." &&
$2 != "./." &&
$2 !~ /^[0-9]+\/[0-9]+$/ {
    print
}' \
> "${PREFIX}.bad_GT_tokens.txt"

echo "==========================================="
echo "CHECKING COLUMN CONSISTENCY"
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

echo "==========================================="
echo "DONE"
echo "==========================================="

echo "Final VCF:"
echo "${PREFIX}.cleaned.vcf.gz"
