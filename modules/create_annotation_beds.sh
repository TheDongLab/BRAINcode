#!/bin/bash
#SBATCH --job-name=make_annotation_beds
#SBATCH --output=make_annotation_beds.out
#SBATCH --error=make_annotation_beds.err
#SBATCH --time=6-00:00:00
#SBATCH -p week
#SBATCH --mem=100G
#SBATCH --cpus-per-task=6

######################
# Date: 2/16/2026
# Version: 1.2
######################

set -euo pipefail
module load BEDTools

ANNOT_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode
STAR_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR
FASTA_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta

GTF=${ANNOT_DIR}/gencode.v49.annotation.gtf
FAI=${FASTA_DIR}/genome.fa.fai
DFAM_HMM=${STAR_DIR}/dfam_all.hmm
HG38_FA=${FASTA_DIR}/genome.fa
OUTPUT_DIR=${STAR_DIR}/RepeatMasker_output

mkdir -p "${STAR_DIR}"

# -----------------------------
# Genome sizes
# -----------------------------
cut -f1,2 "${FAI}" | sort -k1,1V > "${STAR_DIR}/genome.sizes"

# ----------------------------
# Genes BED6 with geneID__symbol__biotype
# ----------------------------
awk 'BEGIN{OFS="\t"}
    $3=="gene"{
        gene_id=""; symbol=""; biotype="";
        match($0,/gene_id "([^"]+)"/,gid); gene_id=gid[1]
        match($0,/gene_name "([^"]+)"/,gn); symbol=gn[1]
        match($0,/gene_type "([^"]+)"/,gt); biotype=gt[1]
        label=gene_id"__"symbol"__"biotype
        print $1,$4-1,$5,label,0,$7
    }' "${GTF}" \
    | sort -k1,1V -k2,2n \
    > "${STAR_DIR}/genes.sorted.bed"

# ----------------------------
# Exons BED6 (unmerged) with exon IDs
# ----------------------------
awk 'BEGIN{OFS="\t"}
    $3=="exon"{
        gene_id=""; symbol=""; biotype=""; exon_id=""
        match($0,/gene_id "([^"]+)"/,gid); gene_id=gid[1]
        match($0,/gene_name "([^"]+)"/,gn); symbol=gn[1]
        match($0,/gene_type "([^"]+)"/,gt); biotype=gt[1]
        match($0,/exon_id "([^"]+)"/,eid); exon_id=eid[1]
        label=gene_id"__"symbol"__"biotype"__"exon_id
        print $1,$4-1,$5,label,0,$7
    }' "${GTF}" \
    | sort -k1,1V -k2,2n -k3,3n -k4,4 -k6,6 \
    | uniq \
    > "${STAR_DIR}/exons.unmerged.bed"

# ----------------------------
# Exons merged per gene
# ----------------------------
sort -k1,1V -k6,6 -k2,2n "${STAR_DIR}/exons.unmerged.bed" \
| awk 'BEGIN{OFS="\t"}
{
    split($4,a,"__")
    gene=a[1]
    key=$1 FS gene FS $6

    if(!(key in current_end)){
        cid[key]=1
        cluster_start[key,1]=$2
        cluster_end[key,1]=$3
    }

    if($2 <= cluster_end[key,cid[key]]){
        if($3 > cluster_end[key,cid[key]])
            cluster_end[key,cid[key]]=$3
    } else {
        cid[key]++
        cluster_start[key,cid[key]]=$2
        cluster_end[key,cid[key]]=$3
    }

    # enforce unique labels per cluster
    label_key=key SUBSEP cid[key] SUBSEP $4
    if(!(label_key in seen)){
        seen[label_key]=1
        labels[key,cid[key]] = (labels[key,cid[key]] ?
                                labels[key,cid[key]]","$4 : $4)
    }
}
END{
    for(k in labels){
        split(k,a,SUBSEP)
        split(a[1],b,FS)
        chr=b[1]; gene=b[2]; strand=b[3]
        c=a[2]
        print chr, cluster_start[a[1],c],
              cluster_end[a[1],c],
              labels[k], 0, strand
    }
}' \
| sort -k1,1V -k2,2n -k3,3n -k6,6 \
> "${STAR_DIR}/exons.merged.bed"

# ----------------------------
# Introns (unmerged)
# ----------------------------
bedtools subtract \
    -a "${STAR_DIR}/genes.sorted.bed" \
    -b "${STAR_DIR}/exons.unmerged.bed" \
    > "${STAR_DIR}/introns.unmerged.tmp.bed"

# Add transcript label placeholder (use geneID__transcriptX)
awk 'BEGIN{OFS="\t"}{
    label=$4"_intron"
    print $1,$2,$3,label,0,$6
}' "${STAR_DIR}/introns.unmerged.tmp.bed" \
    | sort -k1,1V -k2,2n -k3,3n -k4,4 -k6,6 \
    | uniq \
    > "${STAR_DIR}/introns.unmerged.bed"

rm -f "${STAR_DIR}/introns.unmerged.tmp.bed"

# ----------------------------
# Introns merged per gene
# ----------------------------
sort -k1,1V -k6,6 -k2,2n "${STAR_DIR}/introns.unmerged.bed" \
| awk 'BEGIN{OFS="\t"}
{
    split($4,a,"__")
    gene=a[1]
    key=$1 FS gene FS $6

    if(!(key in current_end)){
        cid[key]=1
        cluster_start[key,1]=$2
        cluster_end[key,1]=$3
    }

    if($2 <= cluster_end[key,cid[key]]){
        if($3 > cluster_end[key,cid[key]])
            cluster_end[key,cid[key]]=$3
    } else {
        cid[key]++
        cluster_start[key,cid[key]]=$2
        cluster_end[key,cid[key]]=$3
    }

    label_key=key SUBSEP cid[key] SUBSEP $4
    if(!(label_key in seen)){
        seen[label_key]=1
        labels[key,cid[key]] = (labels[key,cid[key]] ?
                                labels[key,cid[key]]","$4 : $4)
    }
}
END{
    for(k in labels){
        split(k,a,SUBSEP)
        split(a[1],b,FS)
        chr=b[1]; gene=b[2]; strand=b[3]
        c=a[2]
        print chr, cluster_start[a[1],c],
              cluster_end[a[1],c],
              labels[k], 0, strand
    }
}' \
| sort -k1,1V -k2,2n -k3,3n -k6,6 \
> "${STAR_DIR}/introns.merged.bed"

# -----------------------------
# Intergenic regions
# -----------------------------
bedtools complement -g "${STAR_DIR}/genome.sizes" -i "${STAR_DIR}/genes.sorted.bed" \
    > "${STAR_DIR}/intergenic.tmp.bed"

bedtools sort -g "${STAR_DIR}/genome.sizes" -i "${STAR_DIR}/intergenic.tmp.bed" \
    > "${STAR_DIR}/intergenic.sorted.bed"

rm -f "${STAR_DIR}/intergenic.tmp.bed"

# -----------------------------
# 5kb flanking regions
# -----------------------------
bedtools slop -i "${STAR_DIR}/genes.sorted.bed" -g "${STAR_DIR}/genome.sizes" -b 5000 \
    > "${STAR_DIR}/genes.5kb_flank.tmp.bed"

bedtools sort -g "${STAR_DIR}/genome.sizes" -i "${STAR_DIR}/genes.5kb_flank.tmp.bed" \
    > "${STAR_DIR}/genes.5kb_flank.bed"

rm -f "${STAR_DIR}/genes.5kb_flank.tmp.bed"

# -----------------------------
# Intergenic not near genes
# -----------------------------
bedtools subtract -a "${STAR_DIR}/intergenic.sorted.bed" -b "${STAR_DIR}/genes.5kb_flank.bed" \
    > "${STAR_DIR}/intergenic_notneargene.tmp.bed"

bedtools sort -g "${STAR_DIR}/genome.sizes" -i "${STAR_DIR}/intergenic_notneargene.tmp.bed" \
    > "${STAR_DIR}/intergenic_notneargene.sorted.bed"

rm -f "${STAR_DIR}/intergenic_notneargene.tmp.bed"

# -----------------------------
# rRNA (all nuclear + mitochondrial)
# -----------------------------
awk '$3=="gene" && /rRNA/ {OFS="\t"; print $1,$4-1,$5,$7,$7,$7}' "${GTF}" \
    | sort -k1,1V -k2,2n > "${STAR_DIR}/rRNA.sorted.bed"

# -----------------------------
# Non-rRNA mitochondrial
# -----------------------------
awk '$1=="chrM" && $3=="gene"{OFS="\t";print $1,$4-1,$5}' "${GTF}" \
    | sort -k1,1V -k2,2n > "${STAR_DIR}/mtRNA.sorted.bed"

bedtools subtract -a "${STAR_DIR}/mtRNA.sorted.bed" -b "${STAR_DIR}/rRNA.sorted.bed" \
    > "${STAR_DIR}/non_rRNA_mt.tmp.bed"

bedtools sort -g "${STAR_DIR}/genome.sizes" -i "${STAR_DIR}/non_rRNA_mt.tmp.bed" \
    > "${STAR_DIR}/non_rRNA_mt.sorted.bed"

rm -f "${STAR_DIR}/non_rRNA_mt.tmp.bed"

# -----------------------------
# OLD LINE and SINE extraction (obsoleted 2-9)
# -----------------------------
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
# Scan genome with Dfam HMMs (nhmmer)
# -----------------------------
# mkdir -p "${OUTPUT_DIR}/dfam_hits"

# for cls in LINE SINE ERV; do
#     /home/zw529/donglab/pipelines/modules/hmmer-3.4/bin/nhmmer \
#         --tblout "${OUTPUT_DIR}/${cls}.tbl" \
#         --noali \
#         --dna \
#         --incE 0.01 \
#         "${DFAM_HMM}" \
#         "${HG38_FA}" \
#         | tee "${OUTPUT_DIR}/${cls}.nhmmer.log"
# done

#     # Convert nhmmer --tblout to BED
#     awk 'BEGIN{OFS="\t"}
#          !/^#/ {
#             start=$8-1; end=$9
#             strand=($9>$8)?"+":"-"
#             print $1, start, end, $4, $7, strand
#          }' "${OUTPUT_DIR}/${cls}.tbl" \
#          | sort -k1,1V -k2,2n \
#          | bedtools merge -i - \
#          > "${STAR_DIR}/${cls}.dfam.sorted.bed"
# done
