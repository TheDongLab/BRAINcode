#!/usr/bin/env python3
# ~/donglab/pipelines/scripts/rnaseq/circ_percent_calculation.py

import pandas as pd
import pysam
import sys

if len(sys.argv) != 3:
    print("Usage: python circ_percent_calculation.py <sample_bam> <circularRNA_known.txt>")
    sys.exit(1)

bam_file = sys.argv[1]
circ_file = sys.argv[2]
out_file = circ_file.replace(".txt", "_circ_percentage.txt")

# define the circRNA_known.txt column names
circ_columns = [
    "chr", "start", "end", "circ_ID", "strand_count", "strand",
    "backsplice_start", "backsplice_end", "read_support", "num_junctions",
    "junction_reads", "unknown1", "unknown2", "type", "gene", "transcript",
    "circ_exon_nums", "circ_coords"
]

circ = pd.read_csv(circ_file, sep="\t", header=None, names=circ_columns, dtype=str)

# open BAM file
bam = pysam.AlignmentFile(bam_file, "rb")

def calc_linear_exons(coord_string):
    linear_count = 0
    if pd.isna(coord_string):
        return 0
    parts = coord_string.split("|")
    for p in parts:
        try:
            chr_part, range_part = p.split(":")
            start, end = map(int, range_part.split("-"))
            for read in bam.fetch(chr_part, start, end):
                linear_count += 1
        except ValueError:
            continue
    return linear_count

circ["circ_exon"] = circ["circ_exon_nums"].apply(lambda x: sum(map(int, x.split(","))) if pd.notna(x) else 0)
circ["linear_exon"] = circ["circ_coords"].apply(calc_linear_exons)
circ["circ_percent"] = circ["circ_exon"] / (circ["circ_exon"] + circ["linear_exon"])

# write output with original columns first
circ.to_csv(out_file, sep="\t", index=False)
print(f"[DONE] circRNA percentage written to {out_file}")
