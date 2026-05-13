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
module load PLINK

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
echo "STEP 1: EXTRACT FEMALE chrX"
echo "==========================================="

bcftools view \
-S "${FEMALES_SRC}" \
-r chrX,X \
"${VCF_IN}" \
-Oz \
-o "${PREFIX}.female.raw.vcf.gz"

tabix -f -p vcf "${PREFIX}.female.raw.vcf.gz"

############################################
# STEP 2: INITIAL VARIANT COUNTS
############################################

echo "==========================================="
echo "STEP 2: INITIAL VARIANT COUNTS"
echo "==========================================="

bcftools view -H "${PREFIX}.female.raw.vcf.gz" | wc -l \
> "${PREFIX}.initial_variant_count.txt"

bcftools stats "${PREFIX}.female.raw.vcf.gz" \
> "${PREFIX}.initial.stats.txt"

############################################
# STEP 3: INITIAL GT SUMMARY
############################################

echo "==========================================="
echo "STEP 3: INITIAL GT SUMMARY"
echo "==========================================="

bcftools query \
-f '[%GT\n]' \
"${PREFIX}.female.raw.vcf.gz" \
| sort \
| uniq -c \
> "${PREFIX}.initial_GT_summary.txt"

############################################
# STEP 4: CHECK ORIGINAL ALT TYPES
############################################

echo "==========================================="
echo "STEP 4: ORIGINAL ALT TYPE CHECK"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\n' \
"${PREFIX}.female.raw.vcf.gz" \
| awk '
$4 ~ /</ || $4 == "*" {
    print > "'"${PREFIX}.original_symbolic_alleles.txt"'"
}

length($3)>1 || length($4)>1 || $4 ~ /,/ {
    print > "'"${PREFIX}.original_indels_multiallelics.txt"'"
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
"${PREFIX}.female.raw.vcf.gz" \
-Ou \
> "${PREFIX}.step1.bcf"

############################################
# STEP 6: SPLIT MULTIALLELICS + LEFT NORMALIZE
############################################

echo "==========================================="
echo "STEP 6: SPLIT MULTIALLELICS"
echo "==========================================="

bcftools norm \
-f "${REF}" \
-m -any \
"${PREFIX}.step1.bcf" \
-Ou \
> "${PREFIX}.step2.bcf"

############################################
# STEP 7: POST-NORM DIAGNOSTICS
############################################

echo "==========================================="
echo "STEP 7: POST-NORM DIAGNOSTICS"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\n' \
"${PREFIX}.step2.bcf" \
| awk '
length($3)>1 || length($4)>1 || $4 ~ /,/ {
    print
}
' \
> "${PREFIX}.post_norm_complex_variants.txt"

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
"${PREFIX}.step2.bcf" \
-Ou \
> "${PREFIX}.step3.bcf"

############################################
# STEP 9: REMOVE VERY LARGE INDELS
############################################

echo "==========================================="
echo "STEP 9: REMOVE LARGE INDELS"
echo "==========================================="

bcftools view \
-i 'strlen(REF)<=50 && strlen(ALT)<=50' \
"${PREFIX}.step3.bcf" \
-Ou \
> "${PREFIX}.step4.bcf"

############################################
# STEP 10: REMOVE MONOMORPHIC SITES
############################################

echo "==========================================="
echo "STEP 10: REMOVE MONOMORPHIC"
echo "==========================================="

bcftools view \
-c 1:minor \
"${PREFIX}.step4.bcf" \
-Ou \
> "${PREFIX}.step5.bcf"

############################################
# STEP 11: FIX MISSING GT TOKENS
############################################

echo "==========================================="
echo "STEP 11: FIX MISSING GT"
echo "==========================================="

bcftools +setGT \
"${PREFIX}.step5.bcf" \
-Ou \
-- \
-t q \
-n "./." \
-i 'GT="."' \
> "${PREFIX}.step6.bcf"

############################################
# STEP 12: REMOVE PHASING TAGS
############################################

echo "==========================================="
echo "STEP 12: REMOVE PHASING"
echo "==========================================="

bcftools annotate \
-x FORMAT/PID,FORMAT/PGT \
"${PREFIX}.step6.bcf" \
-Ou \
> "${PREFIX}.step7.bcf"

############################################
# STEP 13: CONVERT PHASED GT TO UNPHASED
############################################

echo "==========================================="
echo "STEP 13: UNPHASE GT"
echo "==========================================="

bcftools view \
"${PREFIX}.step7.bcf" \
| sed 's/|/\//g' \
| bcftools view \
-Ou \
> "${PREFIX}.step8.bcf"

############################################
# STEP 14: FINAL VCF
############################################

echo "==========================================="
echo "STEP 14: WRITE FINAL VCF"
echo "==========================================="

bcftools view \
"${PREFIX}.step8.bcf" \
-Oz \
-o "${PREFIX}.cleaned.vcf.gz"

tabix -f -p vcf "${PREFIX}.cleaned.vcf.gz"

############################################
# STEP 15: FINAL GT SUMMARY
############################################

echo "==========================================="
echo "STEP 15: FINAL GT SUMMARY"
echo "==========================================="

bcftools query \
-f '[%GT\n]' \
"${PREFIX}.cleaned.vcf.gz" \
| sort \
| uniq -c \
> "${PREFIX}.final_GT_summary.txt"

############################################
# STEP 16: FINAL VARIANT COUNTS
############################################

echo "==========================================="
echo "STEP 16: FINAL VARIANT COUNTS"
echo "==========================================="

bcftools view -H "${PREFIX}.cleaned.vcf.gz" | wc -l \
> "${PREFIX}.final_variant_count.txt"

bcftools stats "${PREFIX}.cleaned.vcf.gz" \
> "${PREFIX}.final.stats.txt"

############################################
# STEP 17: BAD GT CHECK
############################################

echo "==========================================="
echo "STEP 17: BAD GT CHECK"
echo "==========================================="

bcftools query \
-f '[%SAMPLE\t%GT\n]' \
"${PREFIX}.cleaned.vcf.gz" \
| awk '
$2 != "." &&
$2 != "./." &&
$2 !~ /^[0-9]+\/[0-9]+$/ {
    print
}
' \
> "${PREFIX}.bad_GT_tokens.txt"

############################################
# STEP 18: COLUMN CONSISTENCY CHECK
############################################

echo "==========================================="
echo "STEP 18: COLUMN CONSISTENCY CHECK"
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
}
' \
> "${PREFIX}.bad_lines.txt"

############################################
# STEP 19: SYMBOLIC ALLELE CHECK
############################################

echo "==========================================="
echo "STEP 19: SYMBOLIC ALLELE CHECK"
echo "==========================================="

bcftools view \
-H \
"${PREFIX}.cleaned.vcf.gz" \
| awk '
$5 ~ /</ || $5 == "*" {
    print
}
' \
> "${PREFIX}.remaining_symbolic_alleles.txt"

############################################
# STEP 20: LARGE INDEL CHECK
############################################

echo "==========================================="
echo "STEP 20: LARGE INDEL CHECK"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\n' \
"${PREFIX}.cleaned.vcf.gz" \
| awk '
length($3)>50 || length($4)>50 {
    print
}
' \
> "${PREFIX}.remaining_large_indels.txt"

############################################
# STEP 21: MULTIALLELIC CHECK
############################################

echo "==========================================="
echo "STEP 21: MULTIALLELIC CHECK"
echo "==========================================="

bcftools query \
-f '%CHROM\t%POS\t%ALT\n' \
"${PREFIX}.cleaned.vcf.gz" \
| awk '
$3 ~ /,/ {
    print
}
' \
> "${PREFIX}.remaining_multiallelic.txt"

############################################
# STEP 22: PLINK VALIDATION
############################################

echo "==========================================="
echo "STEP 22: PLINK VALIDATION"
echo "==========================================="

plink \
--vcf "${PREFIX}.cleaned.vcf.gz" \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--make-bed \
--out "${PREFIX}.plink_test"

############################################
# STEP 23: FINAL SUMMARY
############################################

echo "==========================================="
echo "DONE"
echo "==========================================="

echo ""
echo "FINAL FILE:"
echo "${PREFIX}.cleaned.vcf.gz"

echo ""
echo "QC FILES:"
echo "${PREFIX}.initial_GT_summary.txt"
echo "${PREFIX}.final_GT_summary.txt"
echo "${PREFIX}.bad_GT_tokens.txt"
echo "${PREFIX}.bad_lines.txt"
echo "${PREFIX}.remaining_symbolic_alleles.txt"
echo "${PREFIX}.remaining_large_indels.txt"
echo "${PREFIX}.remaining_multiallelic.txt"
echo "${PREFIX}.initial.stats.txt"
echo "${PREFIX}.final.stats.txt"

echo ""
echo "PLINK PREFIX:"
echo "${PREFIX}.plink_test"
