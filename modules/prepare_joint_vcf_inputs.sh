#!/bin/bash
#SBATCH --job-name=prep_joint_vcf
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/prep_joint_vcf.out   ### change log paths based on experiment
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/prep_joint_vcf.err   ### change log paths based on experiment
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0

# Set to target_ALS or answer_ALS before running
EXPERIMENT="answer_ALS"

BASE="/home/zw529/donglab/data/$EXPERIMENT"
QTL_DIR="$BASE/QTL"
mkdir -p "$QTL_DIR"

VCF_LIST="$QTL_DIR/vcf_merge_list.txt"
SAMPLE_MAP="$QTL_DIR/sample_map.txt"
FINAL_VCF="$QTL_DIR/joint_genotyped_GQ.vcf.gz"

# --- STEP 1: SELECT VCFs AND CREATE GT:GQ SIDECARS ---

if [[ "$EXPERIMENT" == "target_ALS" ]]; then
    echo "Ranking TargetALS VCFs and creating standardized sidecars..."

    find "$BASE" -name "*.vcf.gz" -printf '%h\n' | sort -u | while read -r dir; do
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

        if [[ -n "$BEST_VCF" ]]; then
            SIDECAR="$dir/genotypes_quality.vcf.gz"
            if [[ ! -f "$SIDECAR" ]]; then
                echo "Processing winner in $dir: $(basename "$BEST_VCF")"
                bcftools annotate -x INFO,^FORMAT/GT,FORMAT/GQ \
                    "$BEST_VCF" -O z -o "$SIDECAR"
                tabix -f -p vcf "$SIDECAR"
            fi
        fi
    done

    # --- STEP 2: TARGETALS METADATA RECONCILIATION ---
    python3 - <<'EOF'
import pandas as pd
from pathlib import Path
import subprocess

BASE = Path("/home/zw529/donglab/data/target_ALS")
QTL_DIR = BASE / "QTL"

wgs_meta = pd.read_csv(BASE / "targetALS_wgs_metadata.csv", encoding="utf-8-sig")
rnaseq_meta = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")
wgs_meta.columns = wgs_meta.columns.str.strip().str.lower()
rnaseq_meta.columns = rnaseq_meta.columns.str.strip().str.lower()

shared = set(wgs_meta["externalsubjectid"].dropna()) & \
         set(rnaseq_meta["externalsubjectid"].dropna())
wgs = wgs_meta[wgs_meta["externalsubjectid"].isin(shared)].copy()
sidecars = list(BASE.glob("**/genotypes_quality.vcf.gz"))

def find_sidecar(sample_id):
    sid = str(sample_id).lower()
    sid_alt = sid.replace("-", "_")
    for path in sidecars:
        p = str(path).lower()
        if path.name == "genotypes_quality.vcf.gz" and (sid in p or sid_alt in p):
            return str(path)
    return None

wgs["vcf_path"] = wgs["externalsampleid"].apply(find_sidecar)
found = wgs[wgs["vcf_path"].notna()].drop_duplicates("externalsubjectid")

print(f"Final subject count for merge: {len(found)}")
found["vcf_path"].to_csv(QTL_DIR / "vcf_merge_list.txt", index=False, header=False)

sample_map = []
for _, row in found.iterrows():
    result = subprocess.run(
        ["bcftools", "query", "-l", row["vcf_path"]],
        capture_output=True, text=True, check=True
    )
    vcf_id = result.stdout.strip()
    if vcf_id:
        sample_map.append({"vcf_id": vcf_id, "subject_id": row["externalsubjectid"]})

pd.DataFrame(sample_map).to_csv(
    QTL_DIR / "sample_map.txt", sep="\t", index=False, header=False
)
EOF

elif [[ "$EXPERIMENT" == "answer_ALS" ]]; then
    echo "Preparing AnswerALS VCFs for RNA-seq participants..."

    RNA_BASE="$BASE/iPSC_Motor_Neuron/RNAseq/Processed"
    WGS_BASE="$BASE/Blood/WGS/Processed"
    SIDECAR_DIR="$QTL_DIR/sidecars"
    CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"

    mkdir -p "$SIDECAR_DIR"
    : > "$VCF_LIST"
    : > "$SAMPLE_MAP"

    for rna_dir in "$RNA_BASE"/*/; do
        [[ -d "$rna_dir" ]] || continue
        donor=$(basename "$rna_dir")
        candidates=("$WGS_BASE/$donor"/*recalibrated.haplotypeCalls.vcf.gz)

        if ((${#candidates[@]} != 1)); then
            echo "[WARN] $donor: expected one unannotated recalibrated VCF; found ${#candidates[@]}. Skipping."
            continue
        fi

        sidecar="$SIDECAR_DIR/$donor.genotypes_quality.vcf.gz"
        bcftools view -r "$CHROMS" "${candidates[0]}" -Ou |
            bcftools annotate -x INFO,^FORMAT/GT,FORMAT/GQ -Oz -o "$sidecar"
        tabix -f -p vcf "$sidecar"

        vcf_id=$(bcftools query -l "$sidecar")
        if [[ -z "$vcf_id" ]]; then
            echo "[WARN] $donor: no VCF sample ID found. Skipping."
            continue
        fi

        printf '%s\n' "$sidecar" >> "$VCF_LIST"
        printf '%s\t%s\n' "$vcf_id" "$donor" >> "$SAMPLE_MAP"
    done

    [[ -s "$VCF_LIST" ]] || {
        echo "[ERROR] No AnswerALS VCFs matched RNA-seq participant folders." >&2
        exit 1
    }

else
    echo "[ERROR] EXPERIMENT must be target_ALS or answer_ALS." >&2
    exit 1
fi

# --- STEP 3: MERGE & REHEADER ---

echo "Merging sidecars (preserving GT and GQ)..."
bcftools merge -f PASS -l "$VCF_LIST" --threads 8 -O z -o "$FINAL_VCF"

echo "Reheadering sequencing IDs to subject IDs..."
bcftools reheader -s "$SAMPLE_MAP" -o "${FINAL_VCF}.tmp" "$FINAL_VCF"
mv "${FINAL_VCF}.tmp" "$FINAL_VCF"
tabix -f -p vcf "$FINAL_VCF"

echo "Processing complete. Joint VCF ready at $FINAL_VCF"
bcftools view -H "$FINAL_VCF" | head -n 1 | cut -f 1-12
