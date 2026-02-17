#!/bin/bash
#SBATCH --job-name=dfam_scan
#SBATCH --output=dfam_scan_%A_%a.out
#SBATCH --error=dfam_scan_%A_%a.err
#SBATCH --time=6-23:59:00
#SBATCH -p week
#SBATCH --mem=100G
#SBATCH --cpus-per-task=6
#SBATCH --array=1-3

#######################
# Date: 2/16/2026
# Version: 1.0
#######################

set -euo pipefail
module load BEDTools

STAR_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR
FASTA_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta
DFAM_HMM=${STAR_DIR}/dfam_all.hmm
HG38_FA=${FASTA_DIR}/genome.fa
OUTPUT_DIR=${STAR_DIR}/RepeatMasker_output

mkdir -p "${OUTPUT_DIR}"

# -----------------------------
# Dfam download
# -----------------------------
# # Download and combine Dfam HMMs (~90GB compressed → ~900GB uncompressed)
# echo "Downloading Dfam HMM chunks..."
# BASE="https://dfam.org/releases/current/families"
# mkdir -p "${STAR_DIR}/dfam_tmp"
# cd "${STAR_DIR}/dfam_tmp"
# for i in {1..10}; do
#   wget -q "$BASE/Dfam-$i.hmm.gz"
# done
#
# echo "Combining Dfam HMMs into dfam_all.hmm..."
# zcat Dfam-*.hmm.gz > "${DFAM_HMM}"
#
# echo "Cleaning up temporary files..."
# cd "${STAR_DIR}"
# rm -rf dfam_tmp

# -----------------------------
# Map SLURM array ID to repeat class
# -----------------------------
case $SLURM_ARRAY_TASK_ID in
    1) cls="LINE" ;;
    2) cls="SINE" ;;
    3) cls="ERV" ;;
    *) echo "Invalid array ID"; exit 1 ;;
esac

# -----------------------------
# Optional: split genome into chromosome-sized chunks for faster scanning
# -----------------------------
CHUNK_DIR="${STAR_DIR}/hg38_chunks"
mkdir -p "${CHUNK_DIR}"

# create a FASTA file per chromosome if not exist
for chrom in $(cut -f1 "${FASTA_DIR}/genome.fa.fai"); do
    CHR_FASTA="${CHUNK_DIR}/${chrom}.fa"
    if [ ! -f "$CHR_FASTA" ]; then
        samtools faidx "${HG38_FA}" "$chrom" > "$CHR_FASTA"
    fi
done

# -----------------------------
# Run nhmmer per chromosome in parallel
# -----------------------------
TMP_BED="${OUTPUT_DIR}/${cls}_tmp_beds"
mkdir -p "$TMP_BED"

for chrom_fa in ${CHUNK_DIR}/*.fa; do
    chr=$(basename "$chrom_fa" .fa)
    (
        /home/zw529/donglab/pipelines/modules/hmmer-3.4/bin/nhmmer \
            --tblout "${TMP_BED}/${chr}.tbl" \
            --noali \
            --cpu 2 \
            --dna \
            --incE 0.01 \
            "${DFAM_HMM}" \
            "$chrom_fa"
    ) &
done

wait

# -----------------------------
# Convert nhmmer tblout to BED and merge per chromosome
# -----------------------------
cat ${TMP_BED}/*.tbl | awk 'BEGIN{OFS="\t"}
    !/^#/ {
        start=$8-1; end=$9
        strand=($9>$8)?"+":"-"
        print $1, start, end, $4, $7, strand
    }' \
    | sort -k1,1V -k2,2n \
    | bedtools merge -i - \
    > "${STAR_DIR}/${cls}.sorted.bed"

# cleanup
rm -rf "$TMP_BED"
