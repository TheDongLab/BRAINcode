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
import subprocess

BASE = Path("/home/zw529/donglab/data/target_ALS")

# 1. Load Metadata with BOM handling and lowercase normalization
wgs_meta    = pd.read_csv(BASE/"targetALS_wgs_metadata.csv", encoding='utf-8-sig')
rnaseq_meta = pd.read_csv(BASE/"targetALS_rnaseq_metadata.csv")

# Standardize column names
wgs_meta.columns    = wgs_meta.columns.str.strip().str.lower()
rnaseq_meta.columns = rnaseq_meta.columns.str.strip().str.lower()

# 2. Identify shared subjects
wgs_subjects = set(wgs_meta["externalsubjectid"].dropna().astype(str).str.strip())
rna_subjects = set(rnaseq_meta["externalsubjectid"].dropna().astype(str).str.strip())
shared_subjects = wgs_subjects.intersection(rna_subjects)

print(f"Unique Subjects in WGS: {len(wgs_subjects)}")
print(f"Unique Subjects in RNAseq: {len(rna_subjects)}")
print(f"Subjects with BOTH (Expected Merge Size): {len(shared_subjects)}")

# 3. Subset WGS metadata to shared subjects
wgs_subset = wgs_meta[wgs_meta["externalsubjectid"].isin(shared_subjects)].copy()

# 4. Robust VCF searching for both CGND and TALS
print("Scanning for VCF files (ignoring SV/CNV/Repeats)...")
all_vcfs = list(BASE.glob("**/*.vcf.gz"))

def find_best_vcf(sid):
    sid = str(sid)
    sid_underscore = sid.replace('-', '_')
    blacklist = ['cnv', 'sv', 'manta', 'canvas', 'repeats', 'g.vcf', 'md5', 'tbi']
    
    matches = []
    for v in all_vcfs:
        v_name = v.name.lower()
        # Match both externalsampleid (e.g., NEUYY221HZA) and externalsubjectid (e.g., CGND-HDA-00169)
        if (sid.lower() in v_name or sid_underscore.lower() in v_name):
            if not any(b in v_name for b in blacklist):
                matches.append(v)
    
    if not matches:
        return None
    
    # Priority: 1. Recalibrated + Annotated (CGND), 2. Hard-Filtered (TALS), 3. Any SNP VCF
    for m in matches:
        if "recalibrated" in m.name and "annotated" in m.name:
            return str(m)
    for m in matches:
        if "hard-filtered" in m.name:
            return str(m)
    return str(matches[0])

wgs_subset["vcf_path"] = wgs_subset["externalsampleid"].apply(find_best_vcf)

found = wgs_subset[wgs_subset["vcf_path"].notna()].drop_duplicates(subset=["vcf_path"])
missing = wgs_subset[wgs_subset["vcf_path"].isna()]

print(f"VCFs successfully matched: {len(found)}")
print(f"VCFs missing: {len(missing)}")

# 5. Extract actual sample names from VCF headers
print("Extracting sample names from VCF headers...")
vcf_sample_map = {}  # Maps VCF internal ID to subject ID
for _, row in found.iterrows():
    vcf_path = row['vcf_path']
    subject_id = row['externalsubjectid']
    try:
        result = subprocess.run(['bcftools', 'query', '-l', vcf_path], 
                              capture_output=True, text=True, check=True)
        vcf_sample_names = result.stdout.strip().split('\n')
        for sample_name in vcf_sample_names:
            vcf_sample_map[sample_name] = subject_id
    except Exception as e:
        print(f"Warning: Could not extract sample from {vcf_path}: {e}")

# 6. Export lists and Sample Mapping for Reheadering
out_dir = BASE / "QTL"
out_dir.mkdir(parents=True, exist_ok=True)

# Merge List for bcftools -l
found["vcf_path"].to_csv(out_dir/"vcf_merge_list.txt", index=False, header=False)

# Tracking Table
found[["externalsubjectid","externalsampleid","vcf_path"]].to_csv(out_dir/"wgs_samples_for_vcf_merge.csv", index=False)

# Sample Map for reheadering: [VCF_Internal_ID] [Subject_ID]
# Now uses actual VCF header names extracted above
mapping_df = pd.DataFrame(list(vcf_sample_map.items()), columns=['vcf_id', 'subject_id'])
mapping_df.to_csv(out_dir/"sample_map.txt", sep="\t", index=False, header=False)

if len(missing) > 0:
    print("\n--- Subjects missing SNP VCFs ---")
    print(missing[["externalsubjectid", "externalsampleid"]].to_string(index=False))

print(f"\nSample mapping contains {len(vcf_sample_map)} entries")
EOF

# 6. Indexing and Merging
VCF_LIST=/home/zw529/donglab/data/target_ALS/QTL/vcf_merge_list.txt
JOINT_VCF=/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped.vcf.gz
SAMPLE_MAP=/home/zw529/donglab/data/target_ALS/QTL/sample_map.txt

echo "Re-indexing VCFs to ensure .tbi existence..."
cat $VCF_LIST | xargs -I {} -P 4 tabix -f -p vcf {}

echo "Merging samples..."
SAMPLE_COUNT=$(wc -l < $VCF_LIST)
echo "Processing $SAMPLE_COUNT VCF files..."

# 1. bcftools merge -i - : Ignores complex INFO merging rules (handles CGND vs TALS INFO differences)
# 2. bcftools annotate -x ^FORMAT/GT : Strips everything except Genotypes (normalizes TALS extra fields: AF, F1R2, F2R1, etc.)
# 3. bcftools reheader : Maps sequencing IDs (NEUYY..., NEUGN...) to subject IDs (CGND..., TALS...)
bcftools merge -f PASS \
    -i - \
    -l $VCF_LIST \
    --threads 4 \
    -O v | bcftools annotate -x ^FORMAT/GT -O z -o $JOINT_VCF

echo "Reheadering to map sequencing IDs to Subject IDs..."
# This maps IDs like NEUXT..., NEUGN... to CGND..., TALS... so genotypes match expression matrix
bcftools reheader -s $SAMPLE_MAP -o ${JOINT_VCF}.tmp $JOINT_VCF
mv ${JOINT_VCF}.tmp $JOINT_VCF
tabix -p vcf $JOINT_VCF

echo "Final Sample Count in VCF:"
bcftools query -l $JOINT_VCF | wc -l

echo "Done. Joint VCF is ready at $JOINT_VCF"
