#!/bin/bash
#SBATCH --job-name=cryptEx_mgr
#SBATCH --output=/home/zw529/donglab/data/target_ALS/CryptEx/run_cryptEx_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/CryptEx/run_cryptEx_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G

set -euo pipefail

# 1. Paths and Global Settings
CRYPTEX_BIN="/home/zw529/donglab/pipelines/modules/rnaseq/bin/CryptEx"
CRYPTEX_SCRIPT="${CRYPTEX_BIN}/cryptex.sh"
MASTER_META="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"

TISSUES=("Motor_Cortex" "Cervical_Spinal_Cord" "Lumbar_Spinal_Cord" "Thoracic_Spinal_Cord" "Frontal_Cortex" "Cerebellum")
SPECIES="human"
PROTEIN="TDP43" 

EXON_GFF_BASE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation"
STAR_REF_DIR="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR"

mkdir -p /home/zw529/donglab/pipelines/scripts/rnaseq/logs

echo "========================================================"
echo "Starting CryptEx Orchestration on TargetALS Cohorts"
echo "Started at: $(date)"
echo "========================================================"

for TISSUE in "${TISSUES[@]}"; do
    echo "--------------------------------------------------------"
    echo "Processing Tissue: ${TISSUE}"
    echo "--------------------------------------------------------"
    
    TISSUE_DIR="/home/zw529/donglab/data/target_ALS/${TISSUE}"
    META_DIR="${TISSUE_DIR}/metadata"
    SUPPORT_FILE="${META_DIR}/${TISSUE}_support.tab"
    
    mkdir -p "$META_DIR"
    
    # 2. Extract entries matching your tissue regex map and use the standardized STAR symlink
    python3 - <<EOF
import pandas as pd
import os

patterns = {
    "Motor_Cortex": "Motor Cortex Lateral|Motor Cortex Medial|Lateral Motor Cortex|Medial Motor Cortex|Primary Motor Cortex L|Primary Motor Cortex M|Cortex_Motor_Unspecified|Cortex_Motor_BA4|BA4 Motor Cortex|Lateral_motor_cortex|Motor Cortex|BA4",
    "Cervical_Spinal_Cord": "Spinal_Cord_Cervical|Cervical Spinal Cord|Cervical_spinal_cord|Spinal_cord_Cervical|Cervical",
    "Lumbar_Spinal_Cord": "Lumbar Spinal Cord|Spinal_Cord_Lumbosacral|Lumbosacral_Spinal_Cord|Lumbar_spinal_cord|Lumbar|Lumbosacral",
    "Thoracic_Spinal_Cord": "Thoracic Spinal Cord|Thoracic",
    "Frontal_Cortex": "Frontal Cortex|Frontal",
    "Cerebellum": "Cerebellum"
}

df = pd.read_csv("${MASTER_META}")

pattern = patterns["${TISSUE}"]
tissue_df = df[df['tissue'].astype(str).str.contains(pattern, case=False, na=False)].copy()

if tissue_df.empty:
    print(f"⚠️ No matches found in master metadata sheet for pattern group: ${TISSUE}")
    exit(0)

support_records = []

for _, row in tissue_df.iterrows():
    raw_sample_id = row['externalsampleid']
    
    # Flip hyphens to underscores exclusively for directory name resolution
    sample_dir = raw_sample_id.replace('-', '_')
    
    # Track the standard STAR output symlink directly
    bam_path = f"/home/zw529/donglab/data/target_ALS/${TISSUE}/RNAseq/Processed/{sample_dir}/STAR.Aligned.sortedByCoord.out.bam"
    
    # Standardize conditions for DEXSeq contrasts
    raw_group = str(row['subject_group'])
    if "Control" in raw_group:
        condition = "Control"
    elif "ALS" in raw_group or "FTLD" in raw_group:
        condition = "Case"
    else:
        condition = "Other"
        
    support_records.append({
        "sample": sample_dir,
        "bam": bam_path,
        "dataset": "${TISSUE}",
        "condition": condition
    })

support_df = pd.DataFrame(support_records)

# Use os.path.exists to verify symlink path target validity before exporting
support_df['file_exists'] = support_df['bam'].apply(os.path.exists)
valid_support = support_df[support_df['file_exists'] == True].drop(columns=['file_exists'])

if valid_support.empty:
    print(f"⚠️ Dropping out: 0 BAM files verified via symlinks for ${TISSUE}.")
else:
    valid_support.to_csv("${SUPPORT_FILE}", sep="\t", index=False)
    print(f"✅ Generated support.tab for ${TISSUE} with {len(valid_support)} active tracks.")
EOF

    # 3. Guard check: Run pipeline block if support file was written successfully
    if [ ! -f "$SUPPORT_FILE" ]; then
        echo "Skipping ${TISSUE} pipeline section due to missing inputs."
        continue
    fi
    
    # 4. Invoke Module Block Execution
    bash "$CRYPTEX_SCRIPT" \
        --species "$SPECIES" \
        --protein "$PROTEIN" \
        --support "$SUPPORT_FILE" \
        --annotation_file "${STAR_REF_DIR}/geneInfo.tab" \
        --gff "${EXON_GFF_BASE}.gtf" \
        --splice_extractor yes \
        --gff_creator yes \
        --read_counter yes \
        --DEXSeq yes \
        --DESeq no \
        --submit no \
        --strict yes \
        --hold_Step1 no
        
    echo "Finished execution block generation for ${TISSUE} setup layout."
done
