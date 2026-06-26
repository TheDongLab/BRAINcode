#!/bin/bash
#SBATCH --job-name=cryptEx_mgr
#SBATCH --output=/home/zw529/donglab/pipelines/scripts/rnaseq/logs/cryptEx_mgr_%j.out
#SBATCH --error=/home/zw529/donglab/pipelines/scripts/rnaseq/logs/cryptEx_mgr_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=4G

set -euo pipefail

# 1. Define Core Module Paths
CRYPTEX_BIN="/home/zw529/donglab/pipelines/modules/rnaseq/bin/CryptEx"
CRYPTEX_SCRIPT="${CRYPTEX_BIN}/cryptex.sh"

# 2. Setup your tissue lists
TISSUES=("Frontal_Cortex" "Motor_Cortex" "Cerebellum" "Cervical_Spinal_Cord" "Lumbar_Spinal_Cord")
SPECIES="human"
PROTEIN="TDP43" # Adjust based on target criteria

# 3. Reference Genome Annotations on Bouchet (Update to actual paths)
EXON_GFF_BASE="/home/zw529/donglab/references/annotations/Homo_sapiens.GRCh38.105_fixed"

echo "========================================================"
echo "Starting CryptEx Orchestration across multiple tissues"
echo "Started at: $(date)"
echo "========================================================"

for TISSUE in "${TISSUES[@]}"; do
    echo "Processing Tissue: ${TISSUE}"
    
    # Define Target paths dynamically based on data structure
    TISSUE_DIR="/home/zw529/donglab/data/target_ALS/${TISSUE}"
    SUPPORT_FILE="${TISSUE_DIR}/metadata/${TISSUE}_support.tab"
    
    # Skip tissue if support file hasn't been generated yet
    if [ ! -f "$SUPPORT_FILE" ]; then
        echo "⚠️ Warning: Support file missing for ${TISSUE} at ${SUPPORT_FILE}. Skipping."
        continue
    fi
    
    # 4. Invoke the underlying script using local paths and bypass its built-in SGE submission
    # We set --submit no because we will wrap individual module tasks directly via Slurm if needed.
    bash "$CRYPTEX_SCRIPT" \
        --species "$SPECIES" \
        --protein "$PROTEIN" \
        --support "$SUPPORT_FILE" \
        --annotation_file "${EXON_GFF_BASE}_annotations.tab" \
        --gff "$EXON_GFF_BASE" \
        --splice_extractor yes \
        --gff_creator yes \
        --read_counter yes \
        --DEXSeq yes \
        --DESeq no \
        --submit no \
        --strict yes \
        --hold_Step1 no
        
    echo "Finished generation blocks for ${TISSUE} layout."
done
