#!/bin/bash
######################
# Date: 2/10/2026
# Version: 1.0
######################

set -euo pipefail

module load BEDTools

ANNOT_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode
STAR_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR
FASTA_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta

GTF=${ANNOT_DIR}/gencode.v49.annotation.gtf
FAI=${FASTA_DIR}/genome.fa.fai
RMSK=${STAR_DIR}/rmsk.txt.gz
DFAM_HMM=${STAR_DIR}/dfam_all.hmm
HG38_FA=${FASTA_DIR}/genome.fa
OUTPUT_DIR=${STAR_DIR}/RepeatMasker_output

mkdir -p "${STAR_DIR}"

# Create genome sizes file
cut -f1,2 "${FAI}" \
    | sort -k1,1V \
    > "${STAR_DIR}/genome.sizes"

# Extract exons
awk '$3=="exon"{OFS="\t";print $1,$4-1,$5}' "${GTF}" \
    | sort -k1,1V -k2,2n \
    > "${STAR_DIR}/exons.sorted.bed"

# Extract genes
awk '$3=="gene"{OFS="\t";print $1,$4-1,$5}' "${GTF}" \
    | sort -k1,1V -k2,2n \
    > "${STAR_DIR}/genes.sorted.bed"

# Create introns
bedtools subtract \
    -a "${STAR_DIR}/genes.sorted.bed" \
    -b "${STAR_DIR}/exons.sorted.bed" \
    > "${STAR_DIR}/introns.tmp.bed"

bedtools sort \
    -g "${STAR_DIR}/genome.sizes" \
    -i "${STAR_DIR}/introns.tmp.bed" \
    > "${STAR_DIR}/introns.sorted.bed"

rm -f "${STAR_DIR}/introns.tmp.bed"

# Copy genes to genic
cp "${STAR_DIR}/genes.sorted.bed" "${STAR_DIR}/genic.sorted.bed"

# Create intergenic regions
bedtools complement \
    -g "${STAR_DIR}/genome.sizes" \
    -i "${STAR_DIR}/genes.sorted.bed" \
    > "${STAR_DIR}/intergenic.tmp.bed"

bedtools sort \
    -g "${STAR_DIR}/genome.sizes" \
    -i "${STAR_DIR}/intergenic.tmp.bed" \
    > "${STAR_DIR}/intergenic.sorted.bed"

rm -f "${STAR_DIR}/intergenic.tmp.bed"

# Create 5kb flanking regions
bedtools slop \
    -i "${STAR_DIR}/genes.sorted.bed" \
    -g "${STAR_DIR}/genome.sizes" \
    -b 5000 \
    > "${STAR_DIR}/genes.5kb_flank.tmp.bed"

bedtools sort \
    -g "${STAR_DIR}/genome.sizes" \
    -i "${STAR_DIR}/genes.5kb_flank.tmp.bed" \
    > "${STAR_DIR}/genes.5kb_flank.bed"

rm -f "${STAR_DIR}/genes.5kb_flank.tmp.bed"

# Intergenic regions not near genes
bedtools subtract \
    -a "${STAR_DIR}/intergenic.sorted.bed" \
    -b "${STAR_DIR}/genes.5kb_flank.bed" \
    > "${STAR_DIR}/intergenic_notneargene.tmp.bed"

bedtools sort \
    -g "${STAR_DIR}/genome.sizes" \
    -i "${STAR_DIR}/intergenic_notneargene.tmp.bed" \
    > "${STAR_DIR}/intergenic_notneargene.sorted.bed"

rm -f "${STAR_DIR}/intergenic_notneargene.tmp.bed"

# Mitochondrial RNA
awk '$1=="chrM" && $3=="gene"{OFS="\t";print $1,$4-1,$5}' "${GTF}" \
    | sort -k1,1V -k2,2n \
    > "${STAR_DIR}/mtRNA.sorted.bed"

# Ribosomal RNA
awk '$1=="chrM" && $0 ~ /rRNA/{OFS="\t";print $1,$4-1,$5}' "${GTF}" \
    | sort -k1,1V -k2,2n \
    > "${STAR_DIR}/rRNA.sorted.bed"

# Non-rRNA mitochondrial
bedtools subtract \
    -a "${STAR_DIR}/mtRNA.sorted.bed" \
    -b "${STAR_DIR}/rRNA.sorted.bed" \
    > "${STAR_DIR}/non_rRNA_mt.tmp.bed"

bedtools sort \
    -g "${STAR_DIR}/genome.sizes" \
    -i "${STAR_DIR}/non_rRNA_mt.tmp.bed" \
    > "${STAR_DIR}/non_rRNA_mt.sorted.bed"

rm -f "${STAR_DIR}/non_rRNA_mt.tmp.bed"

# # Extract LINE elements (repClass contains "LINE")
# zcat "${RMSK}" \
#     | awk 'BEGIN{OFS="\t"} $12 ~ /^LINE/ {print $6, $7, $8}' \
#     | sort -k1,1V -k2,2n \
#     | bedtools merge -i - \
#     > "${STAR_DIR}/LINE.sorted.bed"
# 
# # Extract SINE elements (repClass contains "SINE")
# zcat "${RMSK}" \
#     | awk 'BEGIN{OFS="\t"} $12 ~ /^SINE/ {print $6, $7, $8}' \
#     | sort -k1,1V -k2,2n \
#     | bedtools merge -i - \
#     > "${STAR_DIR}/SINE.sorted.bed"

# Download and combine Dfam HMMs (~90GB compressed → ~900GB uncompressed)
echo "Downloading Dfam HMM chunks..."
BASE="https://dfam.org/releases/current/families"
mkdir -p "${STAR_DIR}/dfam_tmp"
cd "${STAR_DIR}/dfam_tmp"
for i in {1..10}; do
  wget -q "$BASE/Dfam-$i.hmm.gz"
done

echo "Combining Dfam HMMs into dfam_all.hmm..."
zcat Dfam-*.hmm.gz > "${DFAM_HMM}"

echo "Cleaning up temporary files..."
cd "${STAR_DIR}"
rm -rf dfam_tmp

# Run RepeatMasker with Dfam library
echo "Running RepeatMasker with Dfam (this will take 1-2 days)..."
mkdir -p "${OUTPUT_DIR}"
module load RepeatMasker  # Add this line if available on your cluster
RepeatMasker -lib "${DFAM_HMM}" -pa 16 -dir "${OUTPUT_DIR}" "${HG38_FA}"

# New RMSK file from Dfam RepeatMasker run
NEW_RMSK="${OUTPUT_DIR}/$(basename ${HG38_FA}).out"

mkdir -p "${STAR_DIR}"

# Extract LINE elements from Dfam annotations (repClass contains "LINE")
echo "Extracting LINE elements from Dfam annotations..."
zcat "${NEW_RMSK}" \
    | awk 'BEGIN{OFS="\t"} $12 ~ /^LINE/ {print $6, $7, $8}' \
    | sort -k1,1V -k2,2n \
    | bedtools merge -i - \
    > "${STAR_DIR}/LINE.dfam.sorted.bed"

# Extract SINE elements from Dfam annotations (repClass contains "SINE")
echo "Extracting SINE elements from Dfam annotations..."
zcat "${NEW_RMSK}" \
    | awk 'BEGIN{OFS="\t"} $12 ~ /^SINE/ {print $6, $7, $8}' \
    | sort -k1,1V -k2,2n \
    | bedtools merge -i - \
    > "${STAR_DIR}/SINE.dfam.sorted.bed"
