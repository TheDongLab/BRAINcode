#!/bin/bash
#SBATCH --job-name=prepare_joint_vcf
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/prepare_joint_vcf.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/prepare_joint_vcf.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0 Python-bundle-PyPI/2024.06-GCCcore-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0

python3 - <<'EOF'
import pandas as pd
from pathlib import Path

BASE = Path("/home/zw529/donglab/data/target_ALS")

# 1. Load Metadata with BOM handling and lowercase normalization
wgs_meta    = pd.read_csv(BASE/"targetALS_wgs_metadata.csv", encoding='utf-8-sig')
rnaseq_meta = pd.read_csv(BASE/"targetALS_rnaseq_metadata.csv")

# Standardize column names to lowercase and strip whitespace
wgs_meta.columns    = wgs_meta.columns.str.strip().str.lower()
rnaseq_meta.columns = rnaseq_meta.columns.str.strip().str.lower()

# 2. Identify shared subjects
# We use 'externalsubjectid' as the common key
wgs_subjects = set(wgs_meta["externalsubjectid"].dropna().astype(str).str.strip())
rna_subjects = set(rnaseq_meta["externalsubjectid"].dropna().astype(str).str.strip())
shared_subjects = wgs_subjects.intersection(rna_subjects)

print(f"Unique Subjects in WGS: {len(wgs_subjects)}")
print(f"Unique Subjects in RNAseq: {len(rna_subjects)}")
print(f"Subjects with BOTH (Expected Merge Size): {len(shared_subjects)}")

# 3. Subset WGS metadata to shared subjects
wgs_subset = wgs_meta[wgs_meta["externalsubjectid"].isin(shared_subjects)].copy()

# 4. Robust VCF searching
print("Scanning for VCF files (ignoring SV/CNV/Repeats)...")
all_vcfs = list(BASE.glob("**/*.vcf.gz"))

def find_best_vcf(sid):
    sid = str(sid)
    sid_underscore = sid.replace('-', '_')
    
    # Blacklist non-SNP files
    blacklist = ['cnv', 'sv', 'manta', 'canvas', 'repeats', 'g.vcf', 'md5', 'tbi']
    
    matches = []
    for v in all_vcfs:
        v_name = v.name.lower()
        # Check if Sample ID is in filename
        if (sid.lower() in v_name or sid_underscore.lower() in v_name):
            # Ensure it's not a blacklisted structural variant file
            if not any(b in v_name for b in blacklist):
                matches.append(v)
    
    if not matches:
        return None
    
    # Selection Priority: 1. Recalibrated + Annotated, 2. Hard-Filtered, 3. Any SNP VCF
    for m in matches:
        if "recalibrated" in m.name and "annotated" in m.name:
            return str(m)
    for m in matches:
        if "hard-filtered" in m.name:
            return str(m)
    
    return str(matches[0])

# Map the sample IDs to file paths
wgs_subset["vcf_path"] = wgs_subset["externalsampleid"].apply(find_best_vcf)

found = wgs_subset[wgs_subset["vcf_path"].notna()].drop_duplicates(subset=["vcf_path"])
missing = wgs_subset[wgs_subset["vcf_path"].isna()]

print(f"VCFs successfully matched: {len(found)}")
print(f"VCFs missing: {len(missing)}")

# 5. Export lists
out_dir = BASE / "QTL"
out_dir.mkdir(parents=True, exist_ok=True)

found["vcf_path"].to_csv(out_dir/"vcf_merge_list.txt", index=False, header=False)
found[["externalsubjectid","externalsampleid","vcf_path"]].to_csv(out_dir/"wgs_samples_for_vcf_merge.csv", index=False)

if len(missing) > 0:
    print("\n--- Subjects missing SNP VCFs (likely only have SV/CNV) ---")
    print(missing[["externalsubjectid", "externalsampleid"]].to_string(index=False))
EOF

# 6. Indexing and Merging
VCF_LIST=/home/zw529/donglab/data/target_ALS/QTL/vcf_merge_list.txt

echo "Re-indexing new VCFs..."
cat $VCF_LIST | xargs -I {} -P 4 tabix -f -p vcf {}

echo "Merging VCFs into joint_genotyped.vcf.gz..."
bcftools merge -f PASS \
    -l $VCF_LIST \
    --threads 4 \
    -O z \
    -o /home/zw529/donglab/data/target_ALS/QTL/joint_genotyped.vcf.gz

echo "Done. Proceed to PLINK2 QC script."
