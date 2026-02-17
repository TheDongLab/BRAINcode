#!/bin/bash
#SBATCH --job-name=gtfToRefFlat
#SBATCH --output=gtfToRefFlat_%j.out
#SBATCH --error=gtfToRefFlat_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module load Kent_tools

# Convert GTF -> BED12
gtfToGenePred -genePredExt gencode.v49.annotation.gtf gencode.v49.annotation.genepred
genePredToBed gencode.v49.annotation.genepred gencode.v49.annotation.bed12

# Convert BED12 -> refFlat by reordering columns
# BED12 columns: chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
# refFlat columns: geneName, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds
awk 'BEGIN{OFS="\t"} {
    split($11,sizes,","); split($12,starts,",");
    exonStarts=""; exonEnds="";
    for(i=1;i<=$10;i++){
        exonStarts=exonStarts ($2+starts[i-1])",";
        exonEnds=exonEnds ($2+starts[i-1]+sizes[i-1])",";
    }
    print $4,$4,$1,$6,$2,$3,$7,$8,$10,exonStarts,exonEnds
}' gencode.v49.annotation.bed12 > refFlat.comma.txt

# remove trailing commas from blockStarts/blockEnds of refFlat.txt without breaking tabs
awk -F'\t' 'BEGIN{OFS="\t"} {gsub(/,$/,"",$10); gsub(/,$/,"",$11); print}' refFlat.comma.txt > refFlat.txt
