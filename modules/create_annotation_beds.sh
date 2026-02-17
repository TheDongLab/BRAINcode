#!/bin/bash
#SBATCH --job-name=make_annotation_beds
#SBATCH --output=make_annotation_beds.out
#SBATCH --error=make_annotation_beds.err
#SBATCH --time=23:59:00
#SBATCH -p day
#SBATCH --mem=100G
#SBATCH --cpus-per-task=2

######################
# Date: 2/17/2026
# Version: 1.4
######################

set -euo pipefail
module load BEDTools

ANNOT_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode
STAR_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR
FASTA_DIR=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta

GTF=${ANNOT_DIR}/gencode.v49.annotation.gtf
FAI=${FASTA_DIR}/genome.fa.fai
HG38_FA=${FASTA_DIR}/genome.fa

mkdir -p "${STAR_DIR}"

# -----------------------------
# Genome sizes
# -----------------------------
cut -f1,2 "${FAI}" | sort -k1,1V > "${STAR_DIR}/genome.sizes"

# ----------------------------
# Genes BED6
# ----------------------------
awk 'BEGIN{OFS="\t"}
$3=="gene"{
    match($0,/gene_id "([^"]+)"/,gid)
    match($0,/gene_name "([^"]+)"/,gn)
    match($0,/gene_type "([^"]+)"/,gt)
    label=gid[1]"__"gn[1]"__"gt[1]
    print $1,$4-1,$5,label,0,$7
}' "${GTF}" \
| sort -k1,1V -k2,2n -k3,3n -k6,6 \
> "${STAR_DIR}/genes.sorted.bed"

# ----------------------------
# Exons BED6 (unmerged)
# ----------------------------
awk 'BEGIN{OFS="\t"}
$3=="exon"{
    match($0,/gene_id "([^"]+)"/,gid)
    match($0,/gene_name "([^"]+)"/,gn)
    match($0,/gene_type "([^"]+)"/,gt)
    match($0,/exon_id "([^"]+)"/,eid)
    label=gid[1]"__"gn[1]"__"gt[1]"__"eid[1]
    print $1,$4-1,$5,label,0,$7
}' "${GTF}" \
| sort -k1,1V -k2,2n -k3,3n -k6,6 \
| uniq \
> "${STAR_DIR}/exons.unmerged.bed"

# ----------------------------
# Exons merged per gene (correct connectivity)
# ----------------------------

# add gene_id as column 7
awk 'BEGIN{OFS="\t"}{
    split($4,a,"__")
    print $0,a[1]
}' "${STAR_DIR}/exons.unmerged.bed" \
| sort -k1,1V -k6,6 -k7,7 -k2,2n \
| bedtools cluster -s -i - \
| awk 'BEGIN{OFS="\t"}
{
    # columns:
    # 1 chr
    # 2 start
    # 3 end
    # 4 label
    # 5 score
    # 6 strand
    # 7 gene_id
    # 8 cluster_id

    key=$1 FS $6 FS $7 FS $8

    if(!(key in min) || $2 < min[key]) min[key]=$2
    if(!(key in max) || $3 > max[key]) max[key]=$3

    if(!(key SUBSEP $4 in seen)){
        seen[key SUBSEP $4]=1
        labels[key]=(labels[key]?labels[key]","$4:$4)
    }
}
END{
    for(k in labels){
        split(k,a,FS)
        chr=a[1]
        strand=a[2]
        gene=a[3]
        print chr,min[k],max[k],labels[k],0,strand
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
    -s \
> "${STAR_DIR}/introns.unmerged.tmp.bed"

awk 'BEGIN{OFS="\t"}{
    print $1,$2,$3,$4"_intron",0,$6
}' "${STAR_DIR}/introns.unmerged.tmp.bed" \
| sort -k1,1V -k2,2n -k3,3n -k6,6 \
| uniq \
> "${STAR_DIR}/introns.unmerged.bed"

rm -f "${STAR_DIR}/introns.unmerged.tmp.bed"

# ----------------------------
# Introns merged across genes (retain gene + biotype)
# ----------------------------

sort -k1,1V -k6,6 -k2,2n "${STAR_DIR}/introns.unmerged.bed" \
| bedtools cluster -s -i - \
| awk 'BEGIN{OFS="\t"}
{
    # columns:
    # 1 chr
    # 2 start
    # 3 end
    # 4 geneID__symbol__biotype_intron
    # 5 score
    # 6 strand
    # 7 clusterID

    key=$1 FS $6 FS $7   # chr + strand + cluster

    if(!(key in min) || $2 < min[key]) min[key]=$2
    if(!(key in max) || $3 > max[key]) max[key]=$3

    # preserve full label including biotype
    if(!(key SUBSEP $4 in seen)){
        seen[key SUBSEP $4]=1
        labels[key]=(labels[key] ? labels[key]","$4 : $4)
    }
}
END{
    for(k in labels){
        split(k,a,FS)
        chr=a[1]
        strand=a[2]
        print chr,min[k],max[k],labels[k],0,strand
    }
}' \
| sort -k1,1V -k2,2n -k3,3n -k6,6 \
> "${STAR_DIR}/introns.merged.bed"

# ----------------------------
# Remove ANY intron overlapping ANY exon isoform
# ----------------------------

bedtools subtract \
    -a "${STAR_DIR}/introns.merged.bed" \
    -b "${STAR_DIR}/exons.unmerged.bed" \
    -s \
> "${STAR_DIR}/introns.merged.clean.bed"

mv "${STAR_DIR}/introns.merged.clean.bed" \
   "${STAR_DIR}/introns.merged.bed"
   
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
