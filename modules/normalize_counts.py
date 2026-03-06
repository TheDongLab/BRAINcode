#!/usr/bin/env python3
# normalize_counts.py
# Usage: python3 normalize_counts.py <htseqcount.tab> <GTF> <output_dir> --library-type [PE|SE]

import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("htseq_file",     help="HTSeq count table")
parser.add_argument("gtf_file",       help="GTF annotation file")
parser.add_argument("out_dir",        help="Output directory")
parser.add_argument("--library-type", choices=["PE", "SE"], required=True,
                    help="PE (paired-end) or SE (single-end)")
args = parser.parse_args()

htseq_file   = args.htseq_file
gtf_file     = args.gtf_file
out_dir      = args.out_dir
library_type = args.library_type
norm_label   = "FPKM" if library_type == "PE" else "RPKM"

out_file = os.path.join(out_dir, "normalization.tab")

# ── 1. Parse exon lengths from GTF ───────────────────────────────────────
gene_len = {}
with open(gtf_file) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "exon":
            continue
        exon_len = int(fields[4]) - int(fields[3]) + 1
        for part in fields[8].split(";"):
            part = part.strip()
            if part.startswith("gene_id"):
                gene_id = part.split('"')[1]
                gene_len[gene_id] = gene_len.get(gene_id, 0) + exon_len
                break

# ── 2. Parse HTSeq counts ────────────────────────────────────────────────
counts     = {}
gene_order = []
total      = 0
with open(htseq_file) as fh:
    for line in fh:
        if line.startswith("__") or line.startswith("gene_id"):
            continue
        gene_id, count_str = line.rstrip("\n").split("\t")
        cnt = int(count_str)
        counts[gene_id] = cnt
        gene_order.append(gene_id)
        total += cnt

# ── 3. RPK and TPM scaling factor ────────────────────────────────────────
rpk = {}
for g, cnt in counts.items():
    glen = gene_len.get(g, 0)
    rpk[g] = cnt / (glen / 1e3) if glen > 0 else 0.0

sum_rpk   = sum(rpk.values())
scale_tpm = sum_rpk / 1e6 if sum_rpk > 0 else 1.0

# ── 4. Write combined output table ───────────────────────────────────────
with open(out_file, "w") as fout:
    fout.write(f"gene_id\traw_count\tgene_length_bp\tTPM\t{norm_label}\n")
    for g in gene_order:
        cnt      = counts[g]
        glen     = gene_len.get(g, 0)
        tpm_val  = rpk[g] / scale_tpm
        norm_val = (cnt * 1e9) / (glen * total) if (glen > 0 and total > 0) else 0.0
        fout.write(f"{g}\t{cnt}\t{glen}\t{tpm_val:.6f}\t{norm_val:.6f}\n")

sys.stderr.write(f"  Library type       : {library_type}\n")
sys.stderr.write(f"  Normalization      : TPM + {norm_label}\n")
sys.stderr.write(f"  Output             : {out_file}\n")
sys.stderr.write(f"  Total mapped reads : {total:,}\n")
sys.stderr.write(f"  Genes quantified   : {len(counts):,}\n")
