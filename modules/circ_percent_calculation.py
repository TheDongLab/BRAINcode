#!/usr/bin/env python3
# ~/donglab/pipelines/scripts/rnaseq/circ_percent_calculation.py

import pandas as pd
import pysam
import sys

if len(sys.argv) != 4:
    print("Usage: python circ_percent_calculation.py <sample_bam> <circularRNA_known.txt> <bubble.junction.raw>")
    sys.exit(1)

bam_file = sys.argv[1]
circ_file = sys.argv[2]
raw_junction_file = sys.argv[3]
out_file = circ_file.replace(".txt", "_circ_percentage.txt")

# define the circRNA_known.txt column names
circ_columns = [
    "chr",
    "start",
    "end",
    "circ_ID",
    "strand_count",
    "strand",
    "backsplice_start",
    "backsplice_end",
    "read_support",
    "num_junctions",
    "junction_reads",
    "unknown1",
    "unknown2",
    "type",
    "gene",
    "transcript",
    "circ_exon_nums",
    "circ_coords"
]

circ = pd.read_csv(circ_file, sep="\t", header=None, names=circ_columns, dtype=str)

# Load the ORIGINAL (unsnapped) junction coordinates
raw_junc = pd.read_csv(raw_junction_file, sep="\t", header=None, 
                       names=["chr", "start", "end", "name", "score", "strand"])

# Create a mapping from snapped -> original coordinates
# The circ file has snapped coords, we need to map back to original
junction_map = {}
for _, row in raw_junc.iterrows():
    key = f"{row['chr']}:{row['start']}-{row['end']}"
    junction_map[key] = (row['chr'], int(row['start']), int(row['end']))

# open BAM file
bam = pysam.AlignmentFile(bam_file, "rb")


def calc_linear_exons(coord_string, chrom, start, end):
    """
    Use the ORIGINAL coordinates from the raw junction file
    to count linear reads in the flanking exons
    """
    linear_count = 0
    if pd.isna(coord_string) or coord_string == "None":
        return 0
    
    parts = coord_string.split("|")
    for p in parts:
        if p == "None":
            continue
        try:
            chr_part, range_part = p.split(":")
            start_coord, end_coord = map(int, range_part.split("-"))
            for read in bam.fetch(chr_part, start_coord, end_coord):
                linear_count += 1
        except ValueError:
            continue
    
    return linear_count


circ["circ_exon"] = circ["circ_exon_nums"].apply(
    lambda x: sum(map(int, x.split(","))) if pd.notna(x) else 0
)

circ["linear_exon"] = circ.apply(
    lambda row: calc_linear_exons(row["circ_coords"], row["chr"], row["start"], row["end"]),
    axis=1
)

circ["circ_percent"] = circ.apply(
    lambda row: row["circ_exon"] / (row["circ_exon"] + row["linear_exon"]) 
    if (row["circ_exon"] + row["linear_exon"]) > 0 else 1.0,
    axis=1
)

# write output with original columns first
circ.to_csv(out_file, sep="\t", index=False)

print(f"[DONE] circRNA percentage written to {out_file}")
print(f"[INFO] Total circRNAs: {len(circ)}")
print(f"[INFO] circRNAs with linear support: {(circ['linear_exon'].astype(int) > 0).sum()}")
print(f"[INFO] Mean circ_percent: {circ['circ_percent'].astype(float).mean():.3f}")
