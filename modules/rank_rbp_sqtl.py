#!/usr/bin/env python3
"""Print and save RBP sQTL rankings using observed homozygote mean PSI."""
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

TISSUES = [
    "Cerebellum", "Cervical_Spinal_Cord", "Cortex", "Frontal_Cortex",
    "Lumbar_Spinal_Cord", "Motor_Cortex", "Spinal_Cord",
]


def matrix(path, wanted):
    data = {}
    with path.open() as handle:
        samples = handle.readline().split()[1:]
        if len(samples) != len(set(samples)):
            raise ValueError(f"Duplicate subjects: {path}")
        for line in handle:
            key, _, tail = line.rstrip("\n").partition("\t")
            if key not in wanted:
                continue
            if key in data:
                raise ValueError(f"Duplicate matrix row: {key}")
            data[key] = np.array([
                np.nan if x in ("NA", "NaN", "nan", ".", "") else float(x)
                for x in tail.split("\t")
            ])
            if len(data[key]) != len(samples):
                raise ValueError(f"Wrong row length: {key}")
    return samples, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path.home() / "donglab/data/target_ALS")
    args = parser.parse_args()
    root = args.root
    out = root / "QTL/RBP_sQTL_PSI_ranking"
    source = root / "QTL/RBP_significant_QTL.tsv"
    rbp = pd.read_csv(source, sep="\t", dtype=str)
    rbp = rbp[
        (rbp.QTL_type == "sQTL") & rbp.tissue.isin(TISSUES)
        & (pd.to_numeric(rbp.FDR, errors="raise") <= 0.05)
    ]
    genes = {}
    for key, group in rbp.groupby(["tissue", "trait_id", "snpid"], sort=False):
        genes[key] = ";".join(sorted(set(group.Gene_Name.dropna())))

    sigs, locations, targets = {}, {}, set()
    for tissue in TISSUES:
        q = root / tissue / "sQTL"
        s = pd.read_csv(q / "results" / f"{tissue}_sQTL.FDR0.05.txt",
                        sep="\t", dtype={"junction_id": str, "snpid": str})
        s = s[pd.to_numeric(s.FDR, errors="raise") <= 0.05]
        s = s.drop_duplicates(["junction_id", "snpid"])
        s = s[[ (tissue, j, snp) in genes
                for j, snp in zip(s.junction_id, s.snpid) ]].copy() if len(s) else s
        sigs[tissue] = s
        wanted, loc = set(s.snpid), {}
        with (q / "snp_location.txt").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["snpid"] not in wanted:
                    continue
                coord = (row["chr"].removeprefix("chr"), int(row["pos"]))
                if row["snpid"] in loc and loc[row["snpid"]] != coord:
                    raise ValueError(f"Ambiguous SNP location: {tissue} {row['snpid']}")
                loc[row["snpid"]] = coord
                targets.add(coord)
        locations[tissue] = loc

    alleles = defaultdict(set)
    with (root / "QTL/plink/joint_all_chrs_matrixEQTL.raw").open() as handle:
        for token in handle.readline().split()[6:]:
            m = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([^:]+):([^_]+)_(.+)", token)
            if m:
                chrom, pos, ref, alt, counted = m.groups()
                coord = (chrom, int(pos))
                if coord in targets:
                    alleles[coord].add((ref, alt, counted))

    out.mkdir(parents=True, exist_ok=True)
    all_top, audit = [], []
    metrics = [
        "tissue", "RBP_genes", "variant_chr", "variant_pos", "ref", "alt",
        "plink_counted_allele", "n_ref_ref", "n_het", "n_hom_alt",
        "mean_ref_ref_pct", "mean_het_pct", "mean_hom_alt_pct",
        "delta_alt_minus_ref_pp", "abs_delta_pp",
    ]
    display = [
        "top_rank", "RBP_genes", "junction_id", "snpid", "ref", "alt",
        "n_ref_ref", "n_hom_alt", "mean_ref_ref_pct", "mean_hom_alt_pct",
        "delta_alt_minus_ref_pp", "abs_delta_pp", "FDR",
    ]
    for tissue in TISSUES:
        print(f"\nLoading {tissue}...", flush=True)
        q, s = root / tissue / "sQTL", sigs[tissue]
        records, excluded = [], []
        if len(s):
            subjects, psi = matrix(q / f"splicing_{tissue}.txt", set(s.junction_id))
            gs, gts = matrix(q / f"snp_{tissue}.txt", set(s.snpid))
            if set(subjects) != set(gs):
                raise ValueError(f"PSI/genotype subject sets differ: {tissue}")
            gi = {sample: i for i, sample in enumerate(gs)}
            order = [gi[sample] for sample in subjects]
            for row in s.to_dict("records"):
                j, snp = row["junction_id"], row["snpid"]
                coord = locations[tissue].get(snp)
                aa = alleles.get(coord, set())
                reason = None
                if j not in psi or snp not in gts:
                    reason = "missing_matrix_row"
                elif len(aa) != 1:
                    reason = "missing_or_ambiguous_alleles"
                else:
                    ref, alt, counted = next(iter(aa))
                    if counted not in (ref, alt) or "," in alt or ref == alt:
                        reason = "unsupported_alleles"
                if reason:
                    excluded.append(dict(row, tissue=tissue, reason=reason))
                    continue
                x, g = psi[j], gts[snp][order]
                if np.any(np.isfinite(x) & ((x < 0) | (x > 1))):
                    raise ValueError(f"Non-PSI value: {tissue} {j}")
                if np.any(np.isfinite(g) & ~np.isin(g, [0, 1, 2])):
                    raise ValueError(f"Non-hardcall genotype: {tissue} {snp}")
                altg = 2 - g if counted == ref else g
                r = dict(row, tissue=tissue, RBP_genes=genes[tissue, j, snp],
                         variant_chr="chr" + coord[0], variant_pos=coord[1],
                         ref=ref, alt=alt, plink_counted_allele=counted)
                for value, label in [(0, "ref_ref"), (1, "het"), (2, "hom_alt")]:
                    v = x[(altg == value) & np.isfinite(x)]
                    r["n_" + label] = len(v)
                    r["mean_" + label + "_pct"] = 100 * float(v.mean()) if len(v) else np.nan
                r["delta_alt_minus_ref_pp"] = r["mean_hom_alt_pct"] - r["mean_ref_ref_pct"]
                r["abs_delta_pp"] = abs(r["delta_alt_minus_ref_pp"])
                if r["n_ref_ref"] < 5 or r["n_hom_alt"] < 5:
                    excluded.append(dict(r, reason="fewer_than_5_usable_subjects_in_one_or_both_homozygote_groups"))
                else:
                    records.append(r)
        columns = list(dict.fromkeys(list(s.columns) + metrics))
        ranked = pd.DataFrame(records, columns=columns)
        ranked = ranked.sort_values(
            ["abs_delta_pp", "FDR", "junction_id", "snpid"],
            ascending=[False, True, True, True],
        ).reset_index(drop=True)
        ranked.insert(0, "rank", np.arange(1, len(ranked) + 1))
        top = ranked.loc[ranked.groupby("junction_id", sort=False).cumcount() < 2].head(10).copy()
        top.insert(0, "top_rank", np.arange(1, len(top) + 1))
        ranked.to_csv(out / f"{tissue}.ranked_qualifying.tsv", sep="\t", index=False, na_rep="NA")
        top.to_csv(out / f"{tissue}.top10.tsv", sep="\t", index=False, na_rep="NA")
        pd.DataFrame(excluded, columns=columns + ["reason"]).to_csv(
            out / f"{tissue}.excluded.tsv", sep="\t", index=False, na_rep="NA")
        print(f"\n{tissue}: {len(s):,} unique RBP pairs; {len(ranked):,} qualifying; {len(top)} selected")
        print(top[display].to_string(index=False, float_format=lambda v: f"{v:.6g}")
              if len(top) else "No qualifying pairs.", flush=True)
        all_top.append(top)
        audit.append(dict(tissue=tissue, unique_RBP_pairs=len(s),
                          qualifying_pairs=len(ranked), selected_pairs=len(top),
                          excluded_pairs=len(excluded)))
        if len(s):
            del psi, gts
    pd.concat(all_top, ignore_index=True).to_csv(
        out / "top10_per_tissue.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(audit).to_csv(out / "audit.tsv", sep="\t", index=False)
    print(f"\nSaved results to {out}", flush=True)


if __name__ == "__main__":
    main()
