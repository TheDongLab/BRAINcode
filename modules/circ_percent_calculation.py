#!/usr/bin/env python3
# ~/donglab/pipelines/scripts/rnaseq/circ_percent_calculation.py
#
# FIXES vs previous version:
#   1. circ_reads: parsed from FUSIONJUNC_XXXXX/N in back_spliced_junction.bed
#      (previously summed exon index numbers from circ_exon_nums — WRONG)
#   2. linear_reads: summed from junction_reads col 11 in circularRNA_known.txt
#      (previously fetched ALL reads in large flanking windows via BAM — WRONG)
#   3. BAM file no longer needed as input

import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python circ_percent_calculation.py <circularRNA_known.txt> <back_spliced_junction.bed>")
    sys.exit(1)

circ_file         = sys.argv[1]
raw_junction_file = sys.argv[2]
out_file          = circ_file.replace(".txt", "_circ_percentage.txt")

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

raw_junc["circ_reads"] = raw_junc["name"].apply(
    lambda x: int(x.split("/")[1]) if isinstance(x, str) and "/" in x else 0
)

# Build lookup: "chr:start-end" -> circ_read_count
junc_lookup = {}
for _, row in raw_junc.iterrows():
    key = f"{row['chr']}:{row['start']}-{row['end']}"
    junc_lookup[key] = row["circ_reads"]

circ["circ_exon"] = circ.apply(
    lambda row: junc_lookup.get(f"{row['chr']}:{row['start']}-{row['end']}", 0),
    axis=1
)

# -------------------------------------------------------------------------
# Linear reads from col 11 (junction_reads) — CIRCexplorer2 already counted
# reads spanning the internal splice junctions flanking the circRNA.
# Sum comma-separated values per row.
# e.g. "120,87,185" -> 392
# e.g. "1099"       -> 1099
# -------------------------------------------------------------------------
def sum_junction_reads(s):
    if pd.isna(s) or str(s).strip() in ("", "0"):
        return 0
    try:
        return sum(int(x) for x in str(s).split(","))
    except ValueError:
        return 0

circ["linear_exon"] = circ["junction_reads"].apply(sum_junction_reads)

# -------------------------------------------------------------------------
# circ_percent = circ_reads / (circ_reads + linear_reads)
# If both are 0 (e.g. CDR1as), set to 1.0
# -------------------------------------------------------------------------
circ["circ_percent"] = circ.apply(
    lambda row: (
        int(row["circ_exon"]) / (int(row["circ_exon"]) + row["linear_exon"])
        if (int(row["circ_exon"]) + row["linear_exon"]) > 0
        else 1.0
    ),
    axis=1
)

# -------------------------------------------------------------------------
# Write output
# -------------------------------------------------------------------------
circ.to_csv(out_file, sep="\t", index=False)

print(f"[DONE] circRNA percentage written to {out_file}")
print(f"[INFO] Total circRNAs: {len(circ)}")
print(f"[INFO] circRNAs with circ reads > 0: {(circ['circ_exon'].astype(int) > 0).sum()}")
print(f"[INFO] circRNAs with linear support: {(circ['linear_exon'].astype(int) > 0).sum()}")
print(f"[INFO] Mean circ_percent: {circ['circ_percent'].astype(float).mean():.3f}")
