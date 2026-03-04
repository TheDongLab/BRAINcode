#!/usr/bin/env python3
# ~/donglab/pipelines/scripts/rnaseq/circ_percent_calculation.py
#
# circ_reads:   parsed from circular_RNA/N in col 4 of circularRNA_known.txt
#               (authoritative — this is what CIRCexplorer2 used internally)
# linear_reads: summed from junction_reads col 11 of circularRNA_known.txt
#               (reads spanning internal splice junctions flanking the circRNA)
#               single-exon circs (num_junctions==1) get linear=0 by definition
# BAM file not needed

import pandas as pd
import sys

if len(sys.argv) != 2:
    print("Usage: python circ_percent_calculation.py <circularRNA_known.txt>")
    sys.exit(1)

circ_file = sys.argv[1]
out_file  = circ_file.replace(".txt", "_circ_percentage.txt")

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
# circ_reads from col 4: circular_RNA/N -> N
# This is the authoritative back-splice read count from CIRCexplorer2
# -------------------------------------------------------------------------
circ["circ_exon"] = circ["circ_ID"].apply(
    lambda x: int(x.split("/")[1]) if isinstance(x, str) and "/" in x else 0
)

# -------------------------------------------------------------------------
# linear_reads from col 11: sum comma-separated junction read counts
# single-exon circs (num_junctions==1): col 11 is exon size not read count
# so set linear=0 for these
# -------------------------------------------------------------------------
def sum_junction_reads(row):
    if str(row["num_junctions"]).strip() == "1":
        return 0
    s = row["junction_reads"]
    if pd.isna(s) or str(s).strip() in ("", "0"):
        return 0
    try:
        return sum(int(x) for x in str(s).split(","))
    except ValueError:
        return 0

circ["linear_exon"] = circ.apply(sum_junction_reads, axis=1)

# -------------------------------------------------------------------------
# circ_percent = circ_reads / (circ_reads + linear_reads)
# If both are 0, set to 1.0
# -------------------------------------------------------------------------
circ["circ_percent"] = circ.apply(
    lambda row: (
        row["circ_exon"] / (row["circ_exon"] + row["linear_exon"])
        if (row["circ_exon"] + row["linear_exon"]) > 0
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
print(f"[INFO] circRNAs with circ reads > 0: {(circ['circ_exon'] > 0).sum()}")
print(f"[INFO] circRNAs with linear support: {(circ['linear_exon'] > 0).sum()}")
print(f"[INFO] Mean circ_percent: {circ['circ_percent'].mean():.3f}")
