#!/usr/bin/env python3
import pandas as pd
import sys

if len(sys.argv) != 4:
    print("Usage: python circ_percent_calculation.py <circularRNA_known.txt> <ends_counts_file> <output_percentage.txt>")
    sys.exit(1)

circ_file = sys.argv[1]
counts_file = sys.argv[2]
out_file = sys.argv[3]

# 1. Read boundary coverage counts (original bed format + count in column 6)
counts_df = pd.read_csv(counts_file, sep="\t", header=None, 
                        names=["chrom", "start", "end", "circ_id", "score", "strand", "linear_count"])

# Sum start and end boundary count values per circ_id
linear_summed = counts_df.groupby("circ_id")["linear_count"].sum().to_dict()

# 2. Load the circularRNA annotations (no header in this raw output step)
circ_cols = [
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb",
    "exonCount", "exonSizes", "exonOffsets", "readNumber",
    "circType", "geneName", "isoformName", "index", "flankIntron"
]
circ = pd.read_csv(circ_file, sep="\t", header=None, names=circ_cols, dtype=str)

# Map dynamic index keys
circ["circ_id"] = circ["chrom"] + "_" + circ["start"] + "_" + circ["end"]

# Parse the target circ reads from the standard name format (e.g. circ_id/counts)
circ["circ_reads"] = circ["name"].apply(
    lambda x: int(x.split("/")[1]) if isinstance(x, str) and "/" in x else 0
)

# Align the mapped linear counts
circ["linear_reads"] = circ["circ_id"].map(linear_summed).fillna(0).astype(int)

# 3. Calculate: circ / (circ + linear)
circ["circ_percent"] = circ.apply(
    lambda row: (
        row["circ_reads"] / (row["circ_reads"] + row["linear_reads"])
        if (row["circ_reads"] + row["linear_reads"]) > 0
        else 0.0
    ),
    axis=1
)

# 4. Save calculations without the ID helper column (header will be cleanly added by parent pipeline)
circ.drop(columns=["circ_id"], inplace=True)
circ.to_csv(out_file, sep="\t", header=False, index=False)
