#!/bin/bash
#SBATCH --job-name=prepare_joint_vcf
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/prepare_joint_vcf.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS//eQTLprepare_joint_vcf.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0 Python-bundle-PyPI/2024.06-GCCcore-13.3.0
module load BCFtools/1.21-GCC-13.3.0

python3 - <<'EOF'
import pandas as pd
from pathlib import Path

BASE = Path("/home/zw529/donglab/data/target_ALS")

# ── Load metadata
wgs_meta    = pd.read_csv(BASE/"targetALS_wgs_metadata.csv")
rnaseq_meta = pd.read_csv(BASE/"targetALS_rnaseq_metadata.csv")

wgs_meta.columns    = wgs_meta.columns.str.strip()
rnaseq_meta.columns = rnaseq_meta.columns.str.strip()

# ── Step 1: Identify patients with BOTH WGS + RNA-seq
wgs_patients = set(wgs_meta["Externalsubjectid"].dropna())
rna_patients = set(rnaseq_meta["externalsubjectid"].dropna())
shared_patients = wgs_patients.intersection(rna_patients)

print(f"Patients with BOTH WGS + RNAseq: {len(shared_patients)}")

# ── Step 2: Subset WGS metadata to shared patients
wgs_subset = wgs_meta[wgs_meta["Externalsubjectid"].isin(shared_patients)].copy()

# ── Step 3: Locate recalibrated VCFs robustly (pre-scan all VCFs)
all_vcfs = list(BASE.glob("**/*.recalibrated.haplotypeCalls.vcf.gz"))
wgs_subset["vcf_path"] = wgs_subset["Externalsampleid"].apply(
    lambda sid: next(
        (str(v) for v in all_vcfs if sid.replace('-', '_') in v.name or sid in v.name),
        None
    )
)

found = wgs_subset[wgs_subset["vcf_path"].notna()]
found = found.drop_duplicates(subset=["vcf_path"])
missing = wgs_subset[wgs_subset["vcf_path"].isna()]

print(f"\nVCFs found:   {len(found)}")
print(f"VCFs missing: {len(missing)}")

# ── Step 4: Save list for bcftools merge
vcf_list_file = BASE/eQTL/"vcf_merge_list.txt"
unique_vcfs = found["vcf_path"].unique()
with open(vcf_list_file, "w") as f:
    for path in unique_vcfs:
        f.write(path + "\n")

print(f"\nSaved VCF list: {vcf_list_file}")

# ── Step 5: Save tracking table
found[["Externalsubjectid","Externalsampleid","vcf_path"]].to_csv(
    BASE/eQTL/"wgs_samples_for_vcf_merge.csv", index=False
)

# ── Step 6: Report missing
if len(missing) > 0:
    print("\n--- Missing VCFs ---")
    print(missing[["Externalsampleid","Externalsubjectid"]].to_string(index=False))

EOF

# ── Merge VCFs together:
bcftools merge -f PASS -l /home/zw529/donglab/data/target_ALS/vcf_merge_list.txt -O z -o /home/zw529/donglab/data/target_ALS/joint_genotyped.vcf.gz
