#!/bin/bash
#SBATCH --job-name=split_vcf_by_chr
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/split_vcf.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

set -euo pipefail
export LC_ALL=C

module load BCFtools/1.21

# Set to target_ALS or answer_ALS before running.
EXPERIMENT="answer_ALS"

BASE="/home/zw529/donglab/data/${EXPERIMENT}"
BASE_OUTPUT_DIR="${BASE}/QTL"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/chromosome_joint_vcfs"
INPUT_VCF="${BASE_OUTPUT_DIR}/joint_genotyped_GQ.vcf.gz"

mkdir -p "$OUTPUT_DIR"

if [[ "$EXPERIMENT" == "target_ALS" ]]; then
    PSAM_FILE="${BASE_OUTPUT_DIR}/plink/joint_autosomes_filtered.psam"
    METADATA="${BASE}/targetALS_rnaseq_metadata.csv"
    OUTPUT_PREFIX="target_ALS"
    CHROMOSOMES=({1..22} Y)

    # Keep the existing TargetALS sample selection from the PLINK .psam
    awk 'NR>1 {print $1}' "$PSAM_FILE" > "$OUTPUT_DIR/pass_ids.txt"

    # Keep the existing TargetALS sex lists.
    tr -d '\r' < "$METADATA" |
        awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' |
        sort -u > "$OUTPUT_DIR/full_meta.txt"

    grep -Fwf "$OUTPUT_DIR/pass_ids.txt" "$OUTPUT_DIR/full_meta.txt" |
        awk '$2=="male" {print $1}' > "$OUTPUT_DIR/males.txt"
    grep -Fwf "$OUTPUT_DIR/pass_ids.txt" "$OUTPUT_DIR/full_meta.txt" |
        awk '$2=="female" {print $1}' > "$OUTPUT_DIR/females.txt"

    echo "QC samples: $(wc -l < "$OUTPUT_DIR/pass_ids.txt"), Males: $(wc -l < "$OUTPUT_DIR/males.txt"), Females: $(wc -l < "$OUTPUT_DIR/females.txt")"

elif [[ "$EXPERIMENT" == "answer_ALS" ]]; then
    OUTPUT_PREFIX="answer_ALS"
    CHROMOSOMES=({1..22})

    [[ -s "$INPUT_VCF" ]] || {
        echo "[ERROR] Missing joint AnswerALS VCF: $INPUT_VCF" >&2
        exit 1
    }

    # The joint VCF defines the included AnswerALS participants
    bcftools query -l "$INPUT_VCF" > "$OUTPUT_DIR/pass_ids.txt"

    [[ -s "$OUTPUT_DIR/pass_ids.txt" ]] || {
        echo "[ERROR] No samples found in $INPUT_VCF" >&2
        exit 1
    }

    echo "AnswerALS samples in joint VCF: $(wc -l < "$OUTPUT_DIR/pass_ids.txt")"

else
    echo "[ERROR] EXPERIMENT must be target_ALS or answer_ALS." >&2
    exit 1
fi

[[ -s "$INPUT_VCF" ]] || {
    echo "[ERROR] Missing input VCF: $INPUT_VCF" >&2
    exit 1
}

# Region-based extraction requires an index.
if [[ ! -s "${INPUT_VCF}.tbi" && ! -s "${INPUT_VCF}.csi" ]]; then
    bcftools index -f -t "$INPUT_VCF"
fi

# Extract autosomes and Y for TargetALS; autosomes for AnswerALS.
echo "Processing ${EXPERIMENT} chromosomes..."
for chr in "${CHROMOSOMES[@]}"; do
    bcftools view -r "chr${chr}" -S "$OUTPUT_DIR/pass_ids.txt" \
        --force-samples "$INPUT_VCF" -Oz \
        -o "$OUTPUT_DIR/${OUTPUT_PREFIX}_chr${chr}.vcf.gz"
    bcftools index -f -t "$OUTPUT_DIR/${OUTPUT_PREFIX}_chr${chr}.vcf.gz"
done

# Extract chrX and diploidize haploid genotype calls
echo "Processing chrX with diploid conversion..."
RAW_X="$OUTPUT_DIR/${OUTPUT_PREFIX}_chrX_raw.vcf.gz"
DIPLOID_X="$OUTPUT_DIR/${OUTPUT_PREFIX}_chrX_all_diploid.vcf.gz"

bcftools view -r chrX -S "$OUTPUT_DIR/pass_ids.txt" \
    --force-samples "$INPUT_VCF" -Oz -o "$RAW_X"

bcftools view "$RAW_X" |
    awk -F'\t' 'BEGIN { OFS="\t" }
        /^#/ { print; next }
        {
            nfmt = split($9, fmt, ":")
            gt_index = 0
            for (j = 1; j <= nfmt; j++) {
                if (fmt[j] == "GT") {
                    gt_index = j
                    break
                }
            }

            if (gt_index > 0) {
                for (i = 10; i <= NF; i++) {
                    n = split($i, sample, ":")
                    gt = sample[gt_index]

                    # Convert only single-allele calls; leave diploid calls unchanged.
                    if (gt == "." || gt ~ /^[0-9]+$/) {
                        sample[gt_index] = (gt == "." ? "./." : gt "/" gt)
                        rebuilt = sample[1]
                        for (k = 2; k <= n; k++) {
                            rebuilt = rebuilt ":" sample[k]
                        }
                        $i = rebuilt
                    }
                }
            }
            print
        }' |
    bcftools view -Oz -o "$DIPLOID_X"

bcftools index -f -t "$DIPLOID_X"
rm "$RAW_X"

echo "Script complete. Output files are in $OUTPUT_DIR"
