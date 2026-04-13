#!/bin/bash
#SBATCH --job-name=prepare_joint_vcf
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/prepare_joint_vcf.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/prepare_joint_vcf.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

# Load modules
module load Python/3.12.3-GCCcore-13.3.0 Python-bundle-PyPI/2024.06-GCCcore-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0  # Required for tabix

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
# Note: Ensure the ID intersection is clean by stripping whitespace
wgs_patients = set(wgs_meta["Externalsubjectid"].dropna().str.strip())
rna_patients = set(rnaseq_meta["externalsubjectid"].dropna().str.strip())
shared_patients = wgs_patients.intersection(rna_patients)

print(f"Total Unique Subjects in WGS: {len(wgs_patients)}")
print(f"Total Unique Subjects in RNAseq: {len(rna_patients)}")
print(f"Patients with BOTH WGS + RNAseq: {len(shared_patients)}")

# ── Step 2: Subset WGS metadata to shared patients
wgs_subset = wgs_meta[wgs_meta["Externalsubjectid"].str.strip().isin(shared_patients)].copy()

# ── Step 3: Locate VCFs with fuzzy matching
# Broaden glob to find ALL VCFs; we will filter manually for the best match
print("\nScanning directory for all available VCF files...")
all_vcfs = list(BASE.glob("**/*.vcf.gz"))

def find_vcf_robustly(sid):
    sid = str(sid)
    sid_underscore = sid.replace('-', '_')
    # Rank 1: Exact sid or underscore version in name AND 'recalibrated'
    for v in all_vcfs:
        if (sid in v.name or sid_underscore in v.name) and "recalibrated" in v.name:
            return str(v)
    # Rank 2: Fallback to any VCF containing the sid
    for v in all_vcfs:
        if sid in v.name or sid_underscore in v.name:
            return str(v)
    return None

wgs_subset["vcf_path"] = wgs_subset["Externalsampleid"].apply(find_vcf_robustly)

found = wgs_subset[wgs_subset["vcf_path"].notna()]
found = found.drop_duplicates(subset=["vcf_path"])
missing = wgs_subset[wgs_subset["vcf_path"].isna()]

print(f"VCFs found:   {len(found)}")
print(f"VCFs missing: {len(missing)}")

# ── Step 4: Save list for bcftools
vcf_list_file = BASE / "QTL" / "vcf_merge_list.txt"
vcf_list_file.parent.mkdir(parents=True, exist_ok=True)

with open(vcf_list_file, "w") as f:
    for path in found["vcf_path"].unique():
        f.write(path + "\n")

# ── Step 5: Save tracking table
found[["Externalsubjectid","Externalsampleid","vcf_path"]].to_csv(
    BASE / "QTL" / "wgs_samples_for_vcf_merge.csv",
    index=False
)

# ── Step 6: Report the missing ones so we can manually investigate
if len(missing) > 0:
    print("\n--- ATTENTION: MISSING SAMPLES ---")
    print(missing[["Externalsampleid","Externalsubjectid"]].to_string(index=False))
EOF

# ── Step 7: Index and Merge
# Ensure all files are indexed (tabix won't re-index if .tbi already exists)
echo "Verifying indices..."
cat /home/zw529/donglab/data/target_ALS/QTL/vcf_merge_list.txt | xargs -I {} -P 4 tabix -f -p vcf {}

echo "Starting BCFtools merge..."
bcftools merge -f PASS \
    -l /home/zw529/donglab/data/target_ALS/QTL/vcf_merge_list.txt \
    --threads 4 \
    -O z \
    -o /home/zw529/donglab/data/target_ALS/QTL/joint_genotyped.vcf.gz

echo "VCF Re-generation complete."
