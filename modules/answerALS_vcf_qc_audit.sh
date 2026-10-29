#!/bin/bash
#SBATCH --job-name=answerALS_vcf_qc_audit
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/answerALS_vcf_qc_audit.log
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/answerALS_vcf_qc_audit.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

module --force purge
module load BCFtools/1.21-GCC-13.3.0

BASE="/home/zw529/donglab/data/answer_ALS"
QTL_DIR="$BASE/QTL"
PLINK_DIR="$QTL_DIR/plink"
VCF="$QTL_DIR/joint_genotyped_GQ.vcf.gz"
RAW_PSAM="$PLINK_DIR/joint_all_chrs_raw.psam"
RAW_PVAR="$PLINK_DIR/joint_all_chrs_raw.pvar"
VARIANT_FILTERED_PVAR="$PLINK_DIR/temp_variant_filtered.pvar"
FINAL_PSAM="$PLINK_DIR/joint_all_chrs_filtered.psam"
FINAL_PVAR="$PLINK_DIR/joint_all_chrs_filtered.pvar"
MATRIX="$PLINK_DIR/joint_all_chrs_matrixEQTL.raw"
FAIL_LIST="$PLINK_DIR/all_fail_samples.txt"

for required in "$VCF" "$RAW_PSAM" "$RAW_PVAR" "$VARIANT_FILTERED_PVAR" \
                "$FINAL_PSAM" "$FINAL_PVAR" "$MATRIX"; do
    [[ -s "$required" ]] || { echo "[ERROR] Missing or empty required file: $required" >&2; exit 1; }
done

count_pvar() { awk 'NR>1 {n++} END {print n+0}' "$1"; }
count_psam() { awk 'END {print (NR>0 ? NR-1 : 0)}' "$1"; }
psam_ids() {
    awk 'NR==1 {for (i=1;i<=NF;i++) {h=$i; sub(/^#/, "", h); if (h=="IID") col=i} if (!col) col=1; next} NF>=col {print $col}' "$1" | sort -u
}

VCF_SAMPLES=$(bcftools query -l "$VCF" | awk 'END {print NR+0}')
RAW_SAMPLES=$(count_psam "$RAW_PSAM")
FINAL_SAMPLES=$(count_psam "$FINAL_PSAM")
MATRIX_SAMPLES=$(awk 'END {print (NR>0 ? NR-1 : 0)}' "$MATRIX")
RAW_VARIANTS=$(count_pvar "$RAW_PVAR")
VARIANT_FILTERED=$(count_pvar "$VARIANT_FILTERED_PVAR")
FINAL_VARIANTS=$(count_pvar "$FINAL_PVAR")

echo "AnswerALS genotype QC audit (non-imputed joint VCF)"
echo "Input VCF: $VCF"
echo ""
printf '%-40s %s\n' "Samples in joint VCF:" "$VCF_SAMPLES" \
                    "Samples after PLINK import:" "$RAW_SAMPLES" \
                    "Samples in final PLINK set:" "$FINAL_SAMPLES" \
                    "Sample rows in final dosage matrix:" "$MATRIX_SAMPLES"
echo ""
printf '%-40s %s\n' "Variants after PLINK import:" "$RAW_VARIANTS" \
                    "Variants after --geno 0.05:" "$VARIANT_FILTERED" \
                    "Variants in final PLINK set:" "$FINAL_VARIANTS"

if [[ -f "$FAIL_LIST" ]]; then
    FAIL_COUNT=$(awk 'NF {n++} END {print n+0}' "$FAIL_LIST")
    echo "Samples flagged by sex/heterozygosity/relatedness checks: $FAIL_COUNT"
else
    echo "Samples flagged by sex/heterozygosity/relatedness checks: fail list missing"
fi

echo ""
echo "Samples present after import but absent from the final PLINK set:"
comm -23 <(psam_ids "$RAW_PSAM") <(psam_ids "$FINAL_PSAM")

if [[ "$FINAL_SAMPLES" != "$MATRIX_SAMPLES" ]]; then
    echo "[WARN] Final .psam count and dosage-matrix row count differ." >&2
fi
