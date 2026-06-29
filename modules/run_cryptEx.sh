#!/bin/bash
#SBATCH --job-name=cryptEx_mgr
#SBATCH --output=/home/zw529/donglab/data/target_ALS/CryptEx/run_cryptEx_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/CryptEx/run_cryptEx_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=50G

set -euo pipefail

module load R
module load STAR
module load BEDTools
module load SAMtools

# Paths and Global Settings
CRYPTEX_SCRIPT="/home/zw529/donglab/pipelines/modules/rnaseq/bin/CryptEx/cryptex.sh"
MASTER_META="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"

TISSUES=("Motor_Cortex" "Cervical_Spinal_Cord" "Lumbar_Spinal_Cord" "Frontal_Cortex" "Cerebellum")
SPECIES="human"
PROTEIN="TDP43" 

EXON_GFF_BASE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation"
STAR_REF_DIR="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR"

DATA_CRYPTEX_BASE="/home/zw529/donglab/data/target_ALS/CryptEx"

echo "========================================================"
echo "Starting CryptEx Orchestration on TargetALS Cohorts"
echo "Started at: $(date)"
echo "========================================================"

for TISSUE in "${TISSUES[@]}"; do
    echo "--------------------------------------------------------"
    echo "Processing Tissue: ${TISSUE}"
    echo "--------------------------------------------------------"
    
    # Everything for this tissue drops directly into this dedicated output workspace
    TISSUE_OUT_DIR="${DATA_CRYPTEX_BASE}/${TISSUE}"
    mkdir -p "$TISSUE_OUT_DIR"
    
    SUPPORT_FILE="${TISSUE_OUT_DIR}/${TISSUE}_support.tab"
    
    # 2. Use Python to scan all cluster tissue subdirectories that fit the pattern
    python3 - <<EOF
import pandas as pd
import os
import re

patterns = {
    "Motor_Cortex": "Motor Cortex Lateral|Motor Cortex Medial|Lateral Motor Cortex|Medial Motor Cortex|Primary Motor Cortex L|Primary Motor Cortex M|Cortex_Motor_Unspecified|Cortex_Motor_BA4|BA4 Motor Cortex|Lateral_motor_cortex|Motor Cortex|BA4",
    "Cervical_Spinal_Cord": "Spinal_Cord_Cervical|Cervical Spinal Cord|Cervical_spinal_cord|Spinal_cord_Cervical|Cervical",
    "Lumbar_Spinal_Cord": "Lumbar Spinal Cord|Spinal_Cord_Lumbosacral|Lumbosacral_Spinal_Cord|Lumbar_spinal_cord|Lumbar|Lumbosacral",
    "Thoracic_Spinal_Cord": "Thoracic Spinal Cord|Thoracic",
    "Frontal_Cortex": "Frontal Cortex|Frontal",
    "Cerebellum": "Cerebellum"
}

base_data_dir = "/home/zw529/donglab/data/target_ALS"
compiled_regex = re.compile(patterns["${TISSUE}"], re.IGNORECASE)

matching_dirs = []
if os.path.exists(base_data_dir):
    for d in os.listdir(base_data_dir):
        normalized_d = d.replace('_', ' ')
        if compiled_regex.search(normalized_d) or compiled_regex.search(d):
            full_path = os.path.join(base_data_dir, d)
            if os.path.isdir(full_path):
                matching_dirs.append(d)

print(f"🔍 System scanning directories matching '${TISSUE}': {matching_dirs}")

df = pd.read_csv("${MASTER_META}")
pattern = patterns["${TISSUE}"]
tissue_df = df[df['tissue'].astype(str).str.contains(pattern, case=False, na=False)].copy()

if tissue_df.empty:
    print(f"⚠️ No matches found in master metadata sheet for pattern group: ${TISSUE}")
    exit(0)

support_records = []

for _, row in tissue_df.iterrows():
    raw_sample_id = row['externalsampleid']
    sample_underscore = raw_sample_id.replace('-', '_')
    
    final_bam = None
    for folder in matching_dirs:
        path_underscore = f"{base_data_dir}/{folder}/RNAseq/Processed/{sample_underscore}/STAR.Aligned.sortedByCoord.out.bam"
        path_hyphen = f"{base_data_dir}/{folder}/RNAseq/Processed/{raw_sample_id}/STAR.Aligned.sortedByCoord.out.bam"
        
        if os.path.exists(path_underscore):
            final_bam = path_underscore
            sample_name = sample_underscore
            break
        elif os.path.exists(path_hyphen):
            final_bam = path_hyphen
            sample_name = raw_sample_id
            break
            
    if not final_bam:
        continue
        
    raw_group = str(row['subject_group'])
    if "Control" in raw_group:
        condition = "Control"
    elif "ALS" in raw_group or "FTLD" in raw_group:
        condition = "Case"
    else:
        condition = "Other"
        
    support_records.append({
        "sample": sample_name,
        "bam": final_bam,
        "dataset": "${TISSUE}",
        "condition": condition
    })

if len(support_records) == 0:
    print(f"⚠️ Dropping out: Found 0 verified BAM files across checked directories for ${TISSUE}.")
else:
    support_df = pd.DataFrame(support_records)
    support_df.to_csv("${SUPPORT_FILE}", sep="\t", index=False)
    print(f"✅ Generated support.tab for ${TISSUE} with {len(support_df)} validated active tracks.")
EOF

    # 3. Guard check: Run pipeline block if support file was written successfully
    if [ ! -f "$SUPPORT_FILE" ]; then
        echo "Skipping ${TISSUE} pipeline section due to missing inputs."
        continue
    fi
    
    # 4. Invoke the cleaned core module script with our native parameters
    echo "Spawning generation commands for ${TISSUE}..."
    bash "$CRYPTEX_SCRIPT" \
        --species "$SPECIES" \
        --protein "$PROTEIN" \
        --support "$SUPPORT_FILE" \
        --annotation_file "${STAR_REF_DIR}/geneInfo.tab" \
        --gff "${EXON_GFF_BASE}.gtf" \
        --outdir "$TISSUE_OUT_DIR" \
        --paired yes \
        --stranded no \
        --splice_extractor yes \
        --gff_creator yes \
        --read_counter yes \
        --DEXSeq yes \
        --DESeq no \
        --submit yes \
        --strict yes \
        --hold_Step1 no
        
    echo "Finished execution block generation for ${TISSUE} setup layout."
done
