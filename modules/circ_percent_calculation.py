#!/usr/bin/env python3
# ~/donglab/pipelines/scripts/rnaseq/circ_percent_calculation.py
#
# FIXES vs previous version:
#   - circ_reads now correctly parsed from FUSIONJUNC_XXXXX/N in back_spliced_junction.bed
#     (previously was summing exon index numbers from circ_exon_nums column — WRONG)
#   - linear_reads now excludes SA-tagged and supplementary reads
#     (previously counted ALL reads in flanking window including chimeric — WRONG)

import pandas as pd
import pysam
import sys

if len(sys.argv) != 4:
    print("Usage: python circ_percent_calculation.py <sample_bam> <circularRNA_known.txt> <back_spliced_junction.bed>")
    sys.exit(1)

bam_file           = sys.argv[1]
circ_file          = sys.argv[2]
raw_junction_file  = sys.argv[3]
out_file           = circ_file.replace(".txt", "_circ_percentage.txt")

# -------------------------------------------------------------------------
# Column names for circularRNA_known.txt
# -------------------------------------------------------------------------
circ_columns = [
    "chr", "start", "end", "circ_ID", "strand_count", "strand",
    "backsplice_start", "backsplice_end", "read_support",
    "num_junctions", "junction_reads", "unknown1", "unknown2",
    "type", "gene", "transcript", "circ_exon_nums", "circ_coords"
]

circ = pd.read_csv(circ_file, sep="\t", header=None, names=circ_columns, dtype=str)

# -------------------------------------------------------------------------
# Load back_spliced_junction.bed and extract TRUE circ read counts
# Column 4 (name) format: FUSIONJUNC_<id>/<read_count>
# -------------------------------------------------------------------------
raw_junc = pd.read_csv(
    raw_junction_file, sep="\t", header=None,
    names=["chr", "start", "end", "name", "score", "strand"]
)

# Parse read count from FUSIONJUNC name field
raw_junc["circ_reads"] = raw_junc["name"].apply(
    lambda x: int(x.split("/")[1]) if isinstance(x, str) and "/" in x else 0
)

# Build lookup: "chr:start-end" -> circ_read_count
junc_lookup = {}
for _, row in raw_junc.iterrows():
    key = f"{row['chr']}:{row['start']}-{row['end']}"
    junc_lookup[key] = row["circ_reads"]

# Map circ read counts onto circ df using chr:start-end key
circ["circ_exon"] = circ.apply(
    lambda row: junc_lookup.get(f"{row['chr']}:{row['start']}-{row['end']}", 0),
    axis=1
)

# -------------------------------------------------------------------------
# Open BAM and count LINEAR reads in flanking exon windows
# Excludes:
#   - supplementary alignments (chimeric read halves)
#   - SA-tagged reads (back-splice supporting reads)
#   - unmapped reads
# -------------------------------------------------------------------------
bam = pysam.AlignmentFile(bam_file, "rb")

def calc_linear_reads(coord_string):
    """
    Count non-chimeric reads in flanking exon windows from circ_coords column.
    Format: chr:start-end|chr:start-end
    Only counts reads that are:
      - not supplementary
      - not SA-tagged (chimeric)
      - not unmapped
    """
    linear_count = 0
    if pd.isna(coord_string) or coord_string == "None":
        return 0

    for p in coord_string.split("|"):
        if p == "None":
            continue
        try:
            chr_part, range_part = p.split(":")
            start_coord, end_coord = map(int, range_part.split("-"))
            for read in bam.fetch(chr_part, start_coord, end_coord):
                if read.is_unmapped:
                    continue
                if read.is_supplementary:
                    continue
                if read.has_tag("SA"):
                    continue
                linear_count += 1
        except (ValueError, KeyError):
            continue

    return linear_count

circ["linear_exon"] = circ["circ_coords"].apply(calc_linear_reads)

# -------------------------------------------------------------------------
# Calculate circ percent
# If both are 0 (e.g. CDR1as with no linear reads), set to 1.0
# -------------------------------------------------------------------------
circ["circ_percent"] = circ.apply(
    lambda row: (
        row["circ_exon"] / (row["circ_exon"] + row["linear_exon"])
        if (row["circ_exon"] + row["linear_exon"]) > 0
        else 1.0
    ),
    axis=1
)

bam.close()

# -------------------------------------------------------------------------
# Write output
# -------------------------------------------------------------------------
circ.to_csv(out_file, sep="\t", index=False)

print(f"[DONE] circRNA percentage written to {out_file}")
print(f"[INFO] Total circRNAs: {len(circ)}")
print(f"[INFO] circRNAs with circ reads > 0: {(circ['circ_exon'].astype(int) > 0).sum()}")
print(f"[INFO] circRNAs with linear support: {(circ['linear_exon'].astype(int) > 0).sum()}")
print(f"[INFO] Mean circ_percent: {circ['circ_percent'].astype(float).mean():.3f}")
