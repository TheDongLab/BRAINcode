#!/bin/bash
#SBATCH --job-name=TargetALS_ChrX_Robust_fix
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/TargetALS_ChrX_Robust_fix.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/TargetALS_ChrX_Robust_fix.err

set -euo pipefail

############################################
# MODULES
############################################

module purge
module load BCFtools/1.21
module load HTSlib/1.21
module load PLINK/1.9b_7.11-x86_64 

############################################
# INPUTS
############################################

VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
REF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"
USER_DATA_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
FEMALES_SRC="${USER_DATA_DIR}/females.txt"
MALES_SRC="${USER_DATA_DIR}/males.txt"
PREFIX="${OUT_DIR}/target_ALS_chrX_TOPMed"

mkdir -p "${OUT_DIR}"

############################################
# PROCESS FEMALES AND MALES
############################################

for SEX in females males; do
  
  if [ "${SEX}" == "females" ]; then
    SAMPLE_SRC="${FEMALES_SRC}"
    SEX_PREFIX="${PREFIX}.female"
  else
    SAMPLE_SRC="${MALES_SRC}"
    SEX_PREFIX="${PREFIX}.male"
  fi

############################################
# STEP 1: EXTRACT ${SEX} chrX
############################################

echo "==========================================="
echo "STEP 1: EXTRACT ${SEX^^} chrX"
echo "==========================================="

bcftools view \
-S "${SAMPLE_SRC}" \
-r chrX,X \
"${VCF_IN}" \
-Oz \
-o "${SEX_PREFIX}.raw.vcf.gz"

tabix -f -p vcf "${SEX_PREFIX}.raw.vcf.gz"

############################################
# STEP 2: INITIAL VARIANT COUNTS
############################################

echo "==========================================="
echo "STEP 2: INITIAL VARIANT COUNTS"
echo "==========================================="

bcftools view -H "${SEX_PREFIX}.raw.vcf.gz" | wc -l \
> "${SEX_PREFIX}.initial_variant_count.txt"

bcftools stats "${SEX_PREFIX}.raw.vcf.gz" \
> "${SEX_PREFIX}.initial.stats.txt"

############################################
# STEP 3: INITIAL GT SUMMARY
############################################

echo "==========================================="
echo "STEP 3: INITIAL GT SUMMARY"
echo "==========================================="

bcftools query \
-f '[%GT\n]' \
"${SEX_PREFIX}.raw.vcf.gz" \
| sort \
| uniq -c \
> "${SEX_PREFIX}.initial_GT_summary.txt"

############################################
# STEP 4: CHECK ORIGINAL ALT TYPES
############################################

echo "==========================================="
echo "STEP 4: ORIGINAL ALT TYPE CHECK"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\n' \
"${SEX_PREFIX}.raw.vcf.gz" \
| awk '
$4 ~ /</ || $4 == "*" {
    print > "'"${SEX_PREFIX}.original_symbolic_alleles.txt"'"
}

length($3)>1 || length($4)>1 || $4 ~ /,/ {
    print > "'"${SEX_PREFIX}.original_indels_multiallelics.txt"'"
}
'

############################################
# STEP 5: REMOVE SYMBOLIC / NON_REF ALLELES
############################################

echo "==========================================="
echo "STEP 5: REMOVE SYMBOLIC/NON_REF"
echo "==========================================="

bcftools view \
-e 'ALT="<NON_REF>" || ALT="*" || ALT~"<"' \
"${SEX_PREFIX}.raw.vcf.gz" \
-Ou \
> "${SEX_PREFIX}.step1.bcf"

############################################
# STEP 6: SPLIT MULTIALLELICS + LEFT NORMALIZE
############################################

echo "==========================================="
echo "STEP 6: SPLIT MULTIALLELICS"
echo "==========================================="

bcftools norm \
-f "${REF}" \
-m -any \
"${SEX_PREFIX}.step1.bcf" \
-Ou \
> "${SEX_PREFIX}.step2.bcf"

############################################
# STEP 7: POST-NORM DIAGNOSTICS
############################################

echo "==========================================="
echo "STEP 7: POST-NORM DIAGNOSTICS"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\n' \
"${SEX_PREFIX}.step2.bcf" \
| awk '
length($3)>1 || length($4)>1 || $4 ~ /,/ {
    print
}
' \
> "${SEX_PREFIX}.post_norm_complex_variants.txt"

############################################
# STEP 8: KEEP ONLY BIALLELIC SNP/INDEL
############################################

echo "==========================================="
echo "STEP 8: KEEP BIALLELIC SNP/INDEL"
echo "==========================================="

bcftools view \
-m2 \
-M2 \
-v snps,indels \
"${SEX_PREFIX}.step2.bcf" \
-Ou \
> "${SEX_PREFIX}.step3.bcf"

############################################
# STEP 9: REMOVE VERY LARGE INDELS
############################################

echo "==========================================="
echo "STEP 9: REMOVE LARGE INDELS"
echo "==========================================="

bcftools view \
-i 'strlen(REF)<=50 && strlen(ALT)<=50' \
"${SEX_PREFIX}.step3.bcf" \
-Ou \
> "${SEX_PREFIX}.step4.bcf"

############################################
# STEP 10: REMOVE MONOMORPHIC SITES
############################################

echo "==========================================="
echo "STEP 10: REMOVE MONOMORPHIC"
echo "==========================================="

bcftools view \
-c 1:minor \
"${SEX_PREFIX}.step4.bcf" \
-Ou \
> "${SEX_PREFIX}.step5.bcf"

############################################
# STEP 11: MALE HAPLOID VALIDATION (nonPAR)
############################################

if [ "${SEX}" == "males" ]; then

echo "==========================================="
echo "STEP 11: MALE HAPLOID VALIDATION (nonPAR)"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t[%SAMPLE\t%GT\n]' \
"${SEX_PREFIX}.step5.bcf" \
| awk '
($1 == "X" || $1 == "chrX") &&
(($2 >= 2781480 && $2 <= 155701382) || ($2 >= 57217416)) {
    if ($3 !~ /^[0-9]$/ && $3 !~ /^\./ && $3 != "./." && $3 != ".") {
        print
    }
}
' \
> "${SEX_PREFIX}.nonPAR_diploid_errors.txt"

fi

############################################
# STEP 12: FIX MISSING GT TOKENS
############################################

echo "==========================================="
echo "STEP 12: FIX MISSING GT"
echo "==========================================="

bcftools +setGT \
"${SEX_PREFIX}.step5.bcf" \
-Ou \
-- \
-t q \
-n "./." \
-i 'GT="."' \
> "${SEX_PREFIX}.step6.bcf"

############################################
# STEP 13: REMOVE PHASING TAGS
############################################

echo "==========================================="
echo "STEP 13: REMOVE PHASING"
echo "==========================================="

bcftools annotate \
-x FORMAT/PID,FORMAT/PGT \
"${SEX_PREFIX}.step6.bcf" \
-Ou \
> "${SEX_PREFIX}.step7.bcf"

############################################
# STEP 14: CONVERT PHASED GT TO UNPHASED
############################################

echo "==========================================="
echo "STEP 14: UNPHASE GT"
echo "==========================================="

bcftools view \
"${SEX_PREFIX}.step7.bcf" \
| sed 's/|/\//g' \
| bcftools view \
-Ou \
> "${SEX_PREFIX}.step8.bcf"

############################################
# STEP 15: FINAL VCF
############################################

echo "==========================================="
echo "STEP 15: WRITE FINAL VCF"
echo "==========================================="

bcftools view \
"${SEX_PREFIX}.step8.bcf" \
-Oz \
-o "${SEX_PREFIX}.cleaned.vcf.gz"

tabix -f -p vcf "${SEX_PREFIX}.cleaned.vcf.gz"

############################################
# STEP 16: FINAL GT SUMMARY
############################################

echo "==========================================="
echo "STEP 16: FINAL GT SUMMARY"
echo "==========================================="

bcftools query \
-f '[%GT\n]' \
"${SEX_PREFIX}.cleaned.vcf.gz" \
| sort \
| uniq -c \
> "${SEX_PREFIX}.final_GT_summary.txt"

############################################
# STEP 17: FINAL VARIANT COUNTS
############################################

echo "==========================================="
echo "STEP 17: FINAL VARIANT COUNTS"
echo "==========================================="

bcftools view -H "${SEX_PREFIX}.cleaned.vcf.gz" | wc -l \
> "${SEX_PREFIX}.final_variant_count.txt"

bcftools stats "${SEX_PREFIX}.cleaned.vcf.gz" \
> "${SEX_PREFIX}.final.stats.txt"

############################################
# STEP 18: BAD GT CHECK
############################################

echo "==========================================="
echo "STEP 18: BAD GT CHECK"
echo "==========================================="

bcftools query \
-f '[%SAMPLE\t%GT\n]' \
"${SEX_PREFIX}.cleaned.vcf.gz" \
| awk '
$2 != "." &&
$2 != "./." &&
$2 !~ /^[0-9]+\/[0-9]+$/ {
    print
}
' \
> "${SEX_PREFIX}.bad_GT_tokens.txt"

############################################
# STEP 19: COLUMN CONSISTENCY CHECK
############################################

echo "==========================================="
echo "STEP 19: COLUMN CONSISTENCY CHECK"
echo "==========================================="

zcat "${SEX_PREFIX}.cleaned.vcf.gz" \
| awk '
BEGIN{FS="\t"}

/^#CHROM/{
    expected=NF
    next
}

!/^#/ && NF!=expected{
    print "BAD_LINE", NR, NF, expected, $1, $2
}
' \
> "${SEX_PREFIX}.bad_lines.txt"

############################################
# STEP 20: SYMBOLIC ALLELE CHECK
############################################

echo "==========================================="
echo "STEP 20: SYMBOLIC ALLELE CHECK"
echo "==========================================="

bcftools view \
-H \
"${SEX_PREFIX}.cleaned.vcf.gz" \
| awk '
$5 ~ /</ || $5 == "*" {
    print
}
' \
> "${SEX_PREFIX}.remaining_symbolic_alleles.txt"

############################################
# STEP 21: LARGE INDEL CHECK
############################################

echo "==========================================="
echo "STEP 21: LARGE INDEL CHECK"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\n' \
"${SEX_PREFIX}.cleaned.vcf.gz" \
| awk '
length($3)>50 || length($4)>50 {
    print
}
' \
> "${SEX_PREFIX}.remaining_large_indels.txt"

############################################
# STEP 22: MULTIALLELIC CHECK
############################################

echo "==========================================="
echo "STEP 22: MULTIALLELIC CHECK"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t%ALT\n' \
"${SEX_PREFIX}.cleaned.vcf.gz" \
| awk '
$3 ~ /,/ {
    print
}
' \
> "${SEX_PREFIX}.remaining_multiallelic.txt"

############################################
# STEP 23: PLINK VALIDATION
############################################

echo "==========================================="
echo "STEP 23: PLINK VALIDATION"
echo "==========================================="

plink \
--vcf "${SEX_PREFIX}.cleaned.vcf.gz" \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--make-bed \
--out "${SEX_PREFIX}.plink_test"

############################################
# STEP 24: FINAL SUMMARY
############################################

echo "==========================================="
echo "DONE: ${SEX^^} chrX"
echo "==========================================="

echo ""
echo "FINAL FILE:"
echo "${SEX_PREFIX}.cleaned.vcf.gz"

echo ""
echo "QC FILES:"
echo "${SEX_PREFIX}.initial_GT_summary.txt"
echo "${SEX_PREFIX}.final_GT_summary.txt"
echo "${SEX_PREFIX}.bad_GT_tokens.txt"
echo "${SEX_PREFIX}.bad_lines.txt"
echo "${SEX_PREFIX}.remaining_symbolic_alleles.txt"
echo "${SEX_PREFIX}.remaining_large_indels.txt"
echo "${SEX_PREFIX}.remaining_multiallelic.txt"
echo "${SEX_PREFIX}.initial.stats.txt"
echo "${SEX_PREFIX}.final.stats.txt"

echo ""
echo "PLINK PREFIX:"
echo "${SEX_PREFIX}.plink_test"

done
