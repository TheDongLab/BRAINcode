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

PAR="chrX:10001-2781479,chrX:155701383-156030895"

mkdir -p "${OUT_DIR}"

process_group () {

    GROUP="$1"
    SAMPLE_LIST="$2"
    PREFIX="$3"

    echo "=============================="
    echo "${GROUP}"
    echo "=============================="

    bcftools view \
        -S "${SAMPLE_LIST}" \
        -r chrX,X \
        "${VCF_IN}" \
        -Oz \
        -o "${PREFIX}.raw.vcf.gz"

    tabix -f -p vcf "${PREFIX}.raw.vcf.gz"

    bcftools view -H "${PREFIX}.raw.vcf.gz" | wc -l > "${PREFIX}.initial_variant_count.txt"
    bcftools stats "${PREFIX}.raw.vcf.gz" > "${PREFIX}.initial.stats.txt"

    bcftools query -f '[%GT\n]' "${PREFIX}.raw.vcf.gz" | sort | uniq -c > "${PREFIX}.initial_GT_summary.txt"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${PREFIX}.raw.vcf.gz" \
    | awk -v out_sym="${PREFIX}.initial_symbolic.txt" -v out_cmp="${PREFIX}.initial_complex.txt" '
    $4 ~ /<|\\*/ {print > out_sym}
    length($3)>1 || length($4)>1 || $4 ~ /,/ {print > out_cmp}
    '

    bcftools view \
        -e 'ALT="<NON_REF>" || ALT="*" || ALT~"<"' \
        "${PREFIX}.raw.vcf.gz" \
        -Ou > "${PREFIX}.step1.bcf"

    bcftools norm \
        -f "${REF}" \
        -m -any \
        "${PREFIX}.step1.bcf" \
        -Ou > "${PREFIX}.step2.bcf"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${PREFIX}.step2.bcf" \
    | awk -v out="${PREFIX}.post_norm_complex.txt" '
    length($3)>1 || length($4)>1 || $4 ~ /,/ {print > out}
    '

    bcftools view \
        -m2 -M2 -v snps,indels \
        "${PREFIX}.step2.bcf" \
        -Ou > "${PREFIX}.step3.bcf"

    bcftools view \
        -i 'strlen(REF)<=50 && strlen(ALT)<=50' \
        "${PREFIX}.step3.bcf" \
        -Ou > "${PREFIX}.step4.bcf"

    bcftools view \
        -c 1:minor \
        "${PREFIX}.step4.bcf" \
        -Ou > "${PREFIX}.step5.bcf"

    bcftools +setGT \
        "${PREFIX}.step5.bcf" \
        -Ou -- \
        -t q -n "./." -i 'GT="."' \
        > "${PREFIX}.step6.bcf"

    bcftools annotate \
        -x FORMAT/PID,FORMAT/PGT,FORMAT/PS \
        "${PREFIX}.step6.bcf" \
        -Ou > "${PREFIX}.step7.bcf"

    bcftools view "${PREFIX}.step7.bcf" \
    | awk -v par="${PAR}" '
    function in_par(p, a, i, b, c, start, end) {
        n = split(par, a, ",")
        for (i=1; i<=n; i++) {
            split(a[i], b, ":")
            split(b[2], c, "-")
            start=c[1]; end=c[2]
            if (p>=start && p<=end) return 1
        }
        return 0
    }

    BEGIN {FS="\t"}

    /^#/ {print; next}

    {
        pos=$2

        if (in_par(pos)) {
            for (i=10; i<=NF; i++) {
                gsub(/\|/, "/", $i)
            }
        } else {
            for (i=10; i<=NF; i++) {
                if ($i ~ /^[0-9][\/|][0-9]/) {
                    split($i, a, /[\/|]/)
                    $i=a[1]
                }
            }
        }
        print
    }' \
    | bcftools view -Ob -o "${PREFIX}.step8.bcf"

    bcftools view "${PREFIX}.step8.bcf" -Oz -o "${PREFIX}.cleaned.vcf.gz"
    tabix -f -p vcf "${PREFIX}.cleaned.vcf.gz"

    bcftools query -f '[%GT\n]' "${PREFIX}.cleaned.vcf.gz" | sort | uniq -c > "${PREFIX}.final_GT_summary.txt"
    bcftools view -H "${PREFIX}.cleaned.vcf.gz" | wc -l > "${PREFIX}.final_variant_count.txt"
    bcftools stats "${PREFIX}.cleaned.vcf.gz" > "${PREFIX}.final.stats.txt"

    bcftools query -f '[%SAMPLE\t%GT\n]' "${PREFIX}.cleaned.vcf.gz" \
    | awk '$2!="." && $2!="./." && $2!~/^[0-9]+\/[0-9]+$/' \
    > "${PREFIX}.bad_GT_tokens.txt"

    zcat "${PREFIX}.cleaned.vcf.gz" | awk '
    BEGIN{FS="\t"}
    /^#CHROM/{n=NF; next}
    !/^#/ && NF!=n{print "BAD",NR,NF,n,$1,$2}
    ' > "${PREFIX}.bad_lines.txt"

    bcftools view -H "${PREFIX}.cleaned.vcf.gz" \
    | awk '$5 ~ /<|\\*/' > "${PREFIX}.remaining_symbolic.txt"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${PREFIX}.cleaned.vcf.gz" \
    | awk 'length($3)>50 || length($4)>50' > "${PREFIX}.remaining_large_indels.txt"

    bcftools query -f '%CHROM\t%POS\t%ALT\n' "${PREFIX}.cleaned.vcf.gz" \
    | awk '$3 ~ /,/' > "${PREFIX}.remaining_multiallelic.txt"

    plink \
        --vcf "${PREFIX}.cleaned.vcf.gz" \
        --double-id \
        --allow-extra-chr \
        --set-missing-var-ids @:# \
        --make-bed \
        --out "${PREFIX}.plink_test"

    echo "${PREFIX}.cleaned.vcf.gz"
}

process_group "FEMALE" "${FEMALES_SRC}" "${OUT_DIR}/target_ALS_chrX_TOPMed_female"
process_group "MALE" "${MALES_SRC}" "${OUT_DIR}/target_ALS_chrX_TOPMed_male"
