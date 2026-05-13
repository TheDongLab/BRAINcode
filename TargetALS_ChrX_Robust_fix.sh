#!/bin/bash
#SBATCH --job-name=TargetALS_ChrX_Robust_fix
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/TargetALS_ChrX_Robust_fix.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/diagnostics/TargetALS_ChrX_Robust_fix.err

set -euo pipefail

module purge
module load BCFtools/1.21
module load HTSlib/1.21
module load PLINK/1.9b_7.11-x86_64

VCF_IN="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
REF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
OUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/diagnostics"
USER_DATA_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"

FEMALES_SRC="${USER_DATA_DIR}/females.txt"
MALES_SRC="${USER_DATA_DIR}/males.txt"

PREFIX_F="${OUT_DIR}/target_ALS_chrX_TOPMed_female"
PREFIX_M="${OUT_DIR}/target_ALS_chrX_TOPMed_male"

mkdir -p "${OUT_DIR}"

############################################
# PAR REGIONS (hg38)
############################################
PAR_REGIONS="chrX:10001-2781479,chrX:155701383-156030895"

############################################
# FUNCTION
############################################
process_chrX () {

GROUP=$1
SAMPLE_LIST=$2
PREFIX=$3

echo "==========================================="
echo "PROCESSING ${GROUP}"
echo "==========================================="

############################################
# STEP 1: EXTRACT chrX
############################################

bcftools view \
-S "${SAMPLE_LIST}" \
-r chrX,X \
"${VCF_IN}" \
-Oz \
-o "${PREFIX}.raw.vcf.gz"

tabix -f -p vcf "${PREFIX}.raw.vcf.gz"

############################################
# STEP 2: INITIAL VARIANT COUNTS
############################################

bcftools view -H "${PREFIX}.raw.vcf.gz" | wc -l > "${PREFIX}.initial_variant_count.txt"
bcftools stats "${PREFIX}.raw.vcf.gz" > "${PREFIX}.initial.stats.txt"

############################################
# STEP 3: INITIAL GT SUMMARY
############################################

bcftools query -f '[%GT\n]' "${PREFIX}.raw.vcf.gz" \
| sort | uniq -c > "${PREFIX}.initial_GT_summary.txt"

############################################
# STEP 4: ORIGINAL ALT TYPE CHECK
############################################

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${PREFIX}.raw.vcf.gz" \
| awk '
$4 ~ /</ || $4 == "*" {print > "'"${PREFIX}.symbolic_initial.txt"'"}
length($3)>1 || length($4)>1 || $4 ~ /,/ {print > "'"${PREFIX}.complex_initial.txt"'"}
'

############################################
# STEP 5: REMOVE SYMBOLIC / NON_REF
############################################

bcftools view \
-e 'ALT="<NON_REF>" || ALT="*" || ALT~"<"' \
"${PREFIX}.raw.vcf.gz" \
-Ou > "${PREFIX}.step1.bcf"

############################################
# STEP 6: SPLIT MULTIALLELICS + NORMALIZE
############################################

bcftools norm \
-f "${REF}" \
-m -any \
"${PREFIX}.step1.bcf" \
-Ou > "${PREFIX}.step2.bcf"

############################################
# STEP 7: POST-NORM DIAGNOSTICS
############################################

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${PREFIX}.step2.bcf" \
| awk '
length($3)>1 || length($4)>1 || $4 ~ /,/ {print}
' > "${PREFIX}.post_norm_complex.txt"

############################################
# STEP 8: KEEP BIALLELIC SNP/INDEL
############################################

bcftools view \
-m2 -M2 -v snps,indels \
"${PREFIX}.step2.bcf" \
-Ou > "${PREFIX}.step3.bcf"

############################################
# STEP 9: REMOVE LARGE INDELS
############################################

bcftools view \
-i 'strlen(REF)<=50 && strlen(ALT)<=50' \
"${PREFIX}.step3.bcf" \
-Ou > "${PREFIX}.step4.bcf"

############################################
# STEP 10: REMOVE MONOMORPHIC
############################################

bcftools view \
-c 1:minor \
"${PREFIX}.step4.bcf" \
-Ou > "${PREFIX}.step5.bcf"

############################################
# STEP 11: FIX MISSING GT
############################################

bcftools +setGT \
"${PREFIX}.step5.bcf" \
-Ou -- \
-t q -n "./." -i 'GT="."' \
> "${PREFIX}.step6.bcf"

############################################
# STEP 12: REMOVE PHASING TAGS
############################################

bcftools annotate \
-x FORMAT/PID,FORMAT/PGT,FORMAT/PS \
"${PREFIX}.step6.bcf" \
-Ou > "${PREFIX}.step7.bcf"

############################################
# STEP 13: FIX CHR X PLOIDY (CRITICAL)
############################################

# females: diploid everywhere
# males: diploid ONLY in PAR, haploid elsewhere

bcftools view "${PREFIX}.step7.bcf" \
| awk -v par="${PAR_REGIONS}" '
BEGIN{
    split(par, a, ",")
}
function in_par(pos){
    for(i in a){
        split(a[i],b,":")
        split(b[2],c,"-")
        start=c[1]; end=c[2]
        if(pos>=start && pos<=end) return 1
    }
    return 0
}
{
    if($0 ~ /^#/){print; next}
    if($1!="chrX"){print; next}

    pos=$2

    if(in_par(pos)==1){
        gsub(/\|/, "/")
        print
    } else {
        # force haploid outside PAR
        for(i=10;i<=NF;i++){
            if($i ~ /^[0-9][\/|][0-9]/){
                split($i,a,"[\/|]")
                $i=a[1]
            }
        }
        print
    }
}
' | bcftools view -O b > "${PREFIX}.step8.bcf"

############################################
# STEP 14: FINAL VCF
############################################

bcftools view "${PREFIX}.step8.bcf" \
-Oz -o "${PREFIX}.cleaned.vcf.gz"

tabix -f -p vcf "${PREFIX}.cleaned.vcf.gz"

############################################
# STEP 15: FINAL GT SUMMARY
############################################

bcftools query -f '[%GT\n]' "${PREFIX}.cleaned.vcf.gz" \
| sort | uniq -c > "${PREFIX}.final_GT_summary.txt"

############################################
# STEP 16: FINAL VARIANT COUNTS
############################################

bcftools view -H "${PREFIX}.cleaned.vcf.gz" | wc -l > "${PREFIX}.final_variant_count.txt"
bcftools stats "${PREFIX}.cleaned.vcf.gz" > "${PREFIX}.final.stats.txt"

############################################
# STEP 17: BAD GT CHECK
############################################

bcftools query -f '[%SAMPLE\t%GT\n]' "${PREFIX}.cleaned.vcf.gz" \
| awk '$2!="." && $2!="./." && $2!~/^[0-9]+\/[0-9]+$/' \
> "${PREFIX}.bad_GT_tokens.txt"

############################################
# STEP 18: COLUMN CHECK
############################################

zcat "${PREFIX}.cleaned.vcf.gz" | awk '
BEGIN{FS="\t"}
/^#CHROM/{n=NF; next}
!/^#/ && NF!=n{print "BAD",NR,NF,n,$1,$2}
' > "${PREFIX}.bad_lines.txt"

############################################
# STEP 19: SYMBOLIC CHECK
############################################

bcftools view -H "${PREFIX}.cleaned.vcf.gz" \
| awk '$5 ~ /</ || $5=="*"' \
> "${PREFIX}.remaining_symbolic.txt"

############################################
# STEP 20: INDEL CHECK
############################################

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${PREFIX}.cleaned.vcf.gz" \
| awk 'length($3)>50 || length($4)>50' \
> "${PREFIX}.remaining_large_indels.txt"

############################################
# STEP 21: MULTIALLELIC CHECK
############################################

bcftools query -f '%CHROM\t%POS\t%ALT\n' "${PREFIX}.cleaned.vcf.gz" \
| awk '$3 ~ /,/' \
> "${PREFIX}.remaining_multiallelic.txt"

############################################
# STEP 22: PLINK VALIDATION (IMPORTANT: PLINK1.9)
############################################

plink \
--vcf "${PREFIX}.cleaned.vcf.gz" \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--make-bed \
--out "${PREFIX}.plink_test"

############################################
# DONE
############################################

echo "${PREFIX}.cleaned.vcf.gz"
