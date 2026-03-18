#!/bin/bash
#SBATCH --job-name=create_gene_bed6
#SBATCH --output=create_gene_bed6.out
#SBATCH --error=create_gene_bed6.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

set -euo pipefail

module load Kent_tools
module load BEDTools

for GTF in gencode.v49.annotation.gtf gencode.v47.annotation.gtf gencode.v32.annotation.gtf; do
    OUT="${GTF/.gtf/.gene.bed6}"
    TMP="tmp.gene_spans.${GTF%.gtf}.tsv"

    echo "--- Processing $GTF -> $OUT ---"

    # 1. Extract transcript spans per gene -> min txStart / max txEnd per gene
    awk 'BEGIN{OFS="\t"} $3=="transcript"{
        match($0,/gene_id "([^"]+)"/,gid)
        match($0,/gene_type "([^"]+)"/,gt)
        match($0,/gene_name "([^"]+)"/,gn)
        print gid[1], $1, $4-1, $5, $7, gt[1], gn[1]
    }' "$GTF" | sort -k1,1 \
    | awk 'BEGIN{OFS="\t"} {
        gid=$1
        if (!(gid in chr)) {
            chr[gid]=$2; gs[gid]=$3; ge[gid]=$4; str[gid]=$5; gt[gid]=$6; gn[gid]=$7
        } else {
            if ($3 < gs[gid]) gs[gid]=$3
            if ($4 > ge[gid]) ge[gid]=$4
        }
    } END {
        for (gid in chr) print gid, chr[gid], gs[gid], ge[gid], str[gid], gt[gid], gn[gid]
    }' > "$TMP"

    # 2. Collapse to one gene per row as BED6
    awk 'BEGIN{OFS="\t"} {
        gid=$1; chrom=$2; gs=$3; ge=$4; strand=$5; gt=$6; gn=$7
        name = gid"___"gt"___"gn
        print chrom, gs, ge, name, 0, strand
    }' "$TMP" \
    | sort -k1,1V -k2,2n \
    > "$OUT"

    rm "$TMP"

    echo "--- Done ---"
    echo "$(wc -l < "$OUT") genes written to $OUT"
    echo ""
    echo "--- Spot-check: DDX11L16 ---"
    grep "DDX11L16" "$OUT" || echo "(DDX11L16 not found in $OUT)"
    echo ""

done
