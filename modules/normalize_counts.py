#!/usr/bin/env python3
# normalize_counts.py
# Usage: python3 normalize_counts.py <htseqcount.tab> <GTF> <output_dir>

import sys
import os

htseq_file = sys.argv[1]
gtf_file   = sys.argv[2]
out_dir    = sys.argv[3]

tpm_out  = os.path.join(out_dir, "normalization.TPM.tab")
rpkm_out = os.path.join(out_dir, "normalization.RPKM.tab")
fpkm_out = os.path.join(out_dir, "normalization.FPKM.tab")

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

# ── 4. Write output tables ───────────────────────────────────────────────
with open(tpm_out, "w") as ft, open(rpkm_out, "w") as fr, open(fpkm_out, "w") as ff:
    ft.write("gene_id\traw_count\tgene_length_bp\tTPM\n")
    fr.write("gene_id\traw_count\tgene_length_bp\tRPKM\n")
    ff.write("gene_id\traw_count\tgene_length_bp\tFPKM\n")
    for g in gene_order:
        cnt      = counts[g]
        glen     = gene_len.get(g, 0)
        tpm_val  = rpk[g] / scale_tpm
        rpkm_val = (cnt * 1e9) / (glen * total) if (glen > 0 and total > 0) else 0.0
        fpkm_val = rpkm_val
        ft.write(f"{g}\t{cnt}\t{glen}\t{tpm_val:.6f}\n")
        fr.write(f"{g}\t{cnt}\t{glen}\t{rpkm_val:.6f}\n")
        ff.write(f"{g}\t{cnt}\t{glen}\t{fpkm_val:.6f}\n")

sys.stderr.write(f"  Total mapped reads : {total:,}\n")
sys.stderr.write(f"  Genes quantified   : {len(counts):,}\n")
