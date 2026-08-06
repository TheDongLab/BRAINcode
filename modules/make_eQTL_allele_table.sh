#!/usr/bin/env bash
#SBATCH --job-name=make_eQTL_alleles
#SBATCH --output=/home/zw529/donglab/data/GCST90027163/GWAS/make_eQTL_alleles.out
#SBATCH --error=/home/zw529/donglab/data/GCST90027163/GWAS/make_eQTL_alleles.err
#SBATCH --time=9:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail
module purge
module load BCFtools

PLINK_DIR="$HOME/donglab/data/target_ALS/QTL/plink"; EQTL_DIR="$HOME/donglab/data/target_ALS/Cervical_Spinal_Cord/eQTL"
RAW="$PLINK_DIR/joint_all_chrs_matrixEQTL.raw"; BIM="$PLINK_DIR/joint_all_chrs_filtered_bed.bim"
LEADS="$EQTL_DIR/results/Cervical_Spinal_Cord_eQTL.lead_snps.txt"; SNP_LOC="$EQTL_DIR/snp_location.txt"
OUT="$EQTL_DIR/qtl_alleles_GRCh38.tsv"

for f in "$RAW" "$BIM" "$LEADS" "$SNP_LOC"; do [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f"; exit 1; }; done

export RAW BIM LEADS SNP_LOC OUT
python3 <<'PY'
import os, re, sys
import pandas as pd

raw, bim, leads, snp_loc, out = [os.environ[x] for x in ["RAW","BIM","LEADS","SNP_LOC","OUT"]]

# Read only the RAW header. Every genotype column is variantID_COUNTEDALLELE.
with open(raw) as f:
    raw_cols = f.readline().strip().split()

variant_cols = [x for x in raw_cols if x not in {"FID","IID","PAT","MAT","SEX","PHENOTYPE"}]
records = []
for col in variant_cols:
    if "_" not in col:
        records.append((col, None, "missing_counted_allele_suffix"))
        continue
    variant_id, counted = col.rsplit("_", 1)
    records.append((variant_id, counted.upper(), "ok"))

raw_meta = pd.DataFrame(records, columns=["raw_variant_id","effect_allele","raw_parse_status"])

# BIM columns: chromosome, variant ID, genetic distance, position, A1, A2.
bim = pd.read_csv(bim, sep=r"\s+", header=None, names=["bim_chr","bim_id","cm","bim_pos","bim_a1","bim_a2"], dtype=str)
bim["bim_a1"] = bim["bim_a1"].str.upper(); bim["bim_a2"] = bim["bim_a2"].str.upper()

# RAW variant IDs should correspond to BIM IDs before your rsID conversion.
dat = raw_meta.merge(bim, left_on="raw_variant_id", right_on="bim_id", how="left", validate="one_to_one")
dat["other_allele"] = dat.apply(
    lambda r: r.bim_a2 if r.effect_allele == r.bim_a1 else
              r.bim_a1 if r.effect_allele == r.bim_a2 else None, axis=1
)
dat["allele_status"] = dat.apply(
    lambda r: "ok" if r.effect_allele in {r.bim_a1, r.bim_a2} else
              "missing_bim_match" if pd.isna(r.bim_id) else "counted_allele_not_in_bim", axis=1
)

# Your generated SNP-location table contains the final rsIDs and GRCh38 coordinates.
loc = pd.read_csv(snp_loc, sep="\t", dtype=str).rename(columns={"snpid":"final_snpid","chr":"final_chr","pos":"final_pos"})
loc["key"] = loc["final_chr"].str.replace("^chr","",regex=True).str.upper() + ":" + loc["final_pos"]

dat["key"] = dat["bim_chr"].str.replace("^chr","",regex=True).str.upper() + ":" + dat["bim_pos"]
dat = dat.merge(loc[["key","final_snpid","final_chr","final_pos"]], on="key", how="left")

lead = pd.read_csv(leads, sep="\t", dtype=str)
lead_ids = set(lead["snpid"])
dat = dat[dat["final_snpid"].isin(lead_ids)].copy()

outdat = dat[["final_snpid","final_chr","final_pos","effect_allele","other_allele","bim_a1","bim_a2","raw_variant_id","allele_status"]]
outdat.columns = ["snpid","chr","pos","effect_allele","other_allele","bim_A1","bim_A2","raw_variant_id","allele_status"]
outdat = outdat.drop_duplicates("snpid")

outdat.to_csv(out, sep="\t", index=False)

print(f"RAW genotype columns: {len(variant_cols):,}")
print(f"Lead eQTL SNPs requested: {len(lead_ids):,}")
print(f"Lead SNPs with allele records: {len(outdat):,}")
print("\nAllele status:")
print(outdat["allele_status"].value_counts(dropna=False).to_string())

missing = sorted(lead_ids - set(outdat["snpid"]))
print(f"\nLead SNPs missing from allele table: {len(missing):,}")
if missing:
    print("First missing SNPs:", ", ".join(missing[:20]))

bad = outdat[outdat["allele_status"] != "ok"]
if len(bad):
    bad_path = out + ".problems.tsv"
    bad.to_csv(bad_path, sep="\t", index=False)
    print(f"WARNING: {len(bad):,} problematic rows written to {bad_path}")
PY

echo
echo "Output: $OUT"
head -n 5 "$OUT"
