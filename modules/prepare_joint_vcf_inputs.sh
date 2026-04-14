#!/bin/bash
#SBATCH --job-name=prep_joint_vcf
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/prep_joint_vcf.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/prep_joint_vcf.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0

BASE="/home/zw529/donglab/data/target_ALS"
QTL_DIR="$BASE/QTL"
mkdir -p "$QTL_DIR"

# --- STEP 1: RANKING & SIDECAR GENERATION ---
echo "Ranking VCFs and creating standardized sidecars..."

# 1. Identify unique directories
find "$BASE" -name "*.vcf.gz" -printf '%h\n' | sort -u | while read -r dir; do
    
    # 2. Select the "Best" SNP VCF based on priority and size
    BEST_VCF=$(find "$dir" -maxdepth 1 -name "*.vcf.gz" \
        ! -path "*/QTL/*" \
        ! -name "*annotated*" \
        ! -name "*cnv*" \
        ! -name "*sv*" \
        ! -name "*repeats*" \
        ! -name "*manta*" \
        ! -name "*canvas*" \
        ! -name "*.g.vcf.gz" \
        ! -name "genotypes_quality.vcf.gz" | \
        xargs -r ls -S | \
        awk '
            /recalibrated/ { print 1, $0; next }
            /hard-filtered/ { print 2, $0; next }
            { print 3, $0 }
        ' | sort -n | head -n 1 | cut -d' ' -f2-)

    # 3. Create the GT:GQ sidecar for the winner
    if [ -n "$BEST_VCF" ]; then
        SIDECAR="$dir/genotypes_quality.vcf.gz"
        if [ ! -f "$SIDECAR" ]; then
            echo "Processing winner in $dir: $(basename $BEST_VCF)"
            # This step solves the AD length merge error and saves GQ
            bcftools annotate -x INFO,^FORMAT/GT,FORMAT/GQ "$BEST_VCF" -O z -o "$SIDECAR"
            tabix -f -p vcf "$SIDECAR"
        fi
    fi
done

# --- STEP 2: METADATA RECONCILIATION (FIXED) ---
python3 - <<'EOF'
import pandas as pd
from pathlib import Path
import subprocess

BASE = Path("/home/zw529/donglab/data/target_ALS")
QTL_DIR = BASE / "QTL"

# Load Metadata
wgs_meta = pd.read_csv(BASE/"targetALS_wgs_metadata.csv", encoding='utf-8-sig')
rnaseq_meta = pd.read_csv(BASE/"targetALS_rnaseq_metadata.csv")
wgs_meta.columns = wgs_meta.columns.str.strip().str.lower()
rnaseq_meta.columns = rnaseq_meta.columns.str.strip().str.lower()

# Identify shared subjects
shared_subjects = set(wgs_meta["externalsubjectid"].dropna()).intersection(set(rnaseq_meta["externalsubjectid"].dropna()))
wgs_subset = wgs_meta[wgs_meta["externalsubjectid"].isin(shared_subjects)].copy()

# Find ONLY our new sidecar files
all_sidecars = list(BASE.glob("**/genotypes_quality.vcf.gz"))

def find_sidecar(sample_id):
    sid = str(sample_id).lower()
    sid_alt = sid.replace('-', '_')
    for s in all_sidecars:
        s_path = str(s).lower()
        # Strictly match our sidecar filename
        if (sid in s_path or sid_alt in s_path) and s.name == "genotypes_quality.vcf.gz":
            return str(s)
    return None

wgs_subset["vcf_path"] = wgs_subset["externalsampleid"].apply(find_sidecar)

# --- THE CRITICAL FIX ---
# 1. Remove rows where no VCF was found
# 2. Ensure each SUBJECT ID only appears ONCE to prevent bcftools merge duplicates
found = wgs_subset[wgs_subset["vcf_path"].notna()].drop_duplicates(subset=["externalsubjectid"])

print(f"Final subject count for merge (Duplicates removed): {len(found)}")

# Save merge list
found["vcf_path"].to_csv(QTL_DIR/"vcf_merge_list.txt", index=False, header=False)

# Build sample map
vcf_sample_map = []
for _, row in found.iterrows():
    result = subprocess.run(['bcftools', 'query', '-l', row['vcf_path']], capture_output=True, text=True)
    v_id = result.stdout.strip()
    if v_id:
        vcf_sample_map.append({'vcf_id': v_id, 'subject_id': row['externalsubjectid']})

pd.DataFrame(vcf_sample_map).to_csv(QTL_DIR/"sample_map.txt", sep="\t", index=False, header=False)

print(f"Final subject count for merge: {len(found)}")
EOF

# --- STEP 3: MERGE & REHEADER ---
VCF_LIST="$QTL_DIR/vcf_merge_list.txt"
SAMPLE_MAP="$QTL_DIR/sample_map.txt"
FINAL_VCF="$QTL_DIR/joint_genotyped_GQ.vcf.gz"

echo "Merging sidecars (preserving GT and GQ)..."
bcftools merge -f PASS -l "$VCF_LIST" --threads 8 -O z -o "$FINAL_VCF"

echo "Reheadering sequencing IDs to external subject IDs..."
bcftools reheader -s "$SAMPLE_MAP" -o "${FINAL_VCF}.tmp" "$FINAL_VCF"
mv "${FINAL_VCF}.tmp" "$FINAL_VCF"
tabix -f -p vcf "$FINAL_VCF"

echo "Processing complete. Joint VCF ready at $FINAL_VCF"
bcftools view -H "$FINAL_VCF" | head -n 1 | cut -f 1-12
