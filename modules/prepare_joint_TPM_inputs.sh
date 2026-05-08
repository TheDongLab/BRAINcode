#!/bin/bash
#SBATCH --job-name=prep_joint_TPM
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/prep_joint_TPM_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/prep_joint_TPM_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G

set -euo pipefail

if [ $# -eq 0 ]; then
    echo "Usage: sbatch prepare_joint_TPM_inputs.sh <Tissue_Name>"
    exit 1
fi

TISSUE=$1
BASE="/home/zw529/donglab/data/target_ALS"
COV_FILE="$BASE/QTL/covariates.tsv"
OUT_DIR="$BASE/$TISSUE/eQTL"
FINAL_TPM="$OUT_DIR/joint_TPM_matrix.tsv"

mkdir -p "$OUT_DIR"

# --- Environment Setup ---
module purge
module load R/4.4.2-gfbf-2024a
module load Python/3.12.3-GCCcore-13.3.0

echo "Processing Tissue: $TISSUE"

# --- STEP 1: FIND FILES (Fuzzy Path Matching) ---
echo "Searching for files in $BASE..."
find "$BASE" -name "normalization.tab" | grep -i "$TISSUE" | grep "Processed" | grep -v "/eQTL/" > "$OUT_DIR/tpm_file_list.txt" || true

NUM_FILES=$(wc -l < "$OUT_DIR/tpm_file_list.txt")
echo "Found $NUM_FILES potential normalization.tab files."

if [ "$NUM_FILES" -eq 0 ]; then
    echo "Error: No files found for tissue $TISSUE"
    exit 1
fi

# --- STEP 2: PYTHON MERGE & FILTER ---
python3 - <<EOF
import pandas as pd
from pathlib import Path
import sys

pd.set_option('future.no_silent_downcasting', True)

list_file = "$OUT_DIR/tpm_file_list.txt"
output_file = "$FINAL_TPM"
covariate_path = "$COV_FILE"
target_tissue = "$TISSUE"

TISSUE_REMAP = {
    'Motor Cortex Lateral': 'Motor_Cortex', 'Motor Cortex Medial': 'Motor_Cortex',
    'Lateral Motor Cortex': 'Motor_Cortex', 'Medial Motor Cortex': 'Motor_Cortex',
    'Primary Motor Cortex L': 'Motor_Cortex', 'Primary Motor Cortex M': 'Motor_Cortex', 
    'Cortex_Motor_Unspecified': 'Motor_Cortex', 'Cortex_Motor_BA4': 'Motor_Cortex', 
    'BA4 Motor Cortex': 'Motor_Cortex', 'Lateral_motor_cortex': 'Motor_Cortex',
    'Frontal Cortex': 'Frontal_Cortex', 'Cerebellum': 'Cerebellum', 
    'Spinal_Cord_Cervical': 'Cervical_Spinal_Cord', 'Cervical Spinal Cord': 'Cervical_Spinal_Cord', 
    'Cervical_spinal_cord': 'Cervical_Spinal_Cord', 'Spinal_cord_Cervical': 'Cervical_Spinal_Cord', 
    'Lumbar Spinal Cord': 'Lumbar_Spinal_Cord', 'Thoracic Spinal Cord': 'Thoracic_Spinal_Cord', 
    'Spinal_Cord_Lumbosacral': 'Lumbar_Spinal_Cord', 'Lumbosacral_Spinal_Cord': 'Lumbar_Spinal_Cord', 
    'Lumbar_spinal_cord': 'Lumbar_Spinal_Cord'
}

# 1. Load Metadata
meta = pd.read_csv(covariate_path, sep='\t')
meta.columns = meta.columns.str.strip()

id_col = 'externalsampleid'
tissue_col = 'tissue'

# 2. Filter metadata for the target tissue
meta['mapped_tissue'] = meta[tissue_col].map(TISSUE_REMAP).fillna(meta[tissue_col])
meta_subset = meta[meta['mapped_tissue'] == target_tissue].copy()

# 3. Create ID Mapping {folder_format: metadata_format}
# This ensures final headers match the DASHED format in covariates.tsv
id_map = {str(val).replace('-', '_'): str(val) for val in meta_subset[id_col]}
valid_metadata_ids = set(meta_subset[id_col].astype(str))

print(f"Metadata filtered for {target_tissue}: {len(meta_subset)} potential samples.")

# 4. Read file list
try:
    with open(list_file, 'r') as f:
        files = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    sys.exit(1)

all_dfs = []
for fpath in files:
    folder_name = Path(fpath).parent.name
    
    # Resolve the ID to match the metadata format (Dashes)
    meta_id = None
    if folder_name in id_map:
        meta_id = id_map[folder_name]
    elif folder_name.replace('_', '-') in valid_metadata_ids:
        meta_id = folder_name.replace('_', '-')
    elif folder_name in valid_metadata_ids:
        meta_id = folder_name
        
    if meta_id:
        try:
            df = pd.read_csv(fpath, sep='\t', usecols=['gene_id', 'TPM'])
            df['gene_id'] = df['gene_id'].astype(str).str.strip().str.split('.').str[0]
            df = df.drop_duplicates(subset=['gene_id'])
            # Set columns to meta_id so R can find them in the covariates file
            df = df.set_index('gene_id').rename(columns={'TPM': meta_id})
            all_dfs.append(df)
        except Exception as e:
            print(f"Error reading {folder_name}: {e}")

if all_dfs:
    print(f"Merging {len(all_dfs)} dataframes...")
    joint_df = pd.concat(all_dfs, axis=1, join='outer').fillna(0).infer_objects(copy=False)
    
    # Filter: > 0.1 TPM in >= 20% of samples
    threshold = len(all_dfs) * 0.2
    mask = (joint_df > 0.1).sum(axis=1) >= threshold
    joint_df = joint_df[mask]
    
    joint_df.to_csv(output_file, sep='\t')
    print(f"Successfully merged {len(all_dfs)} samples.")
    print(f"Final Matrix dimensions: {joint_df.shape}")
else:
    print("No samples matched metadata IDs.")
    sys.exit(1)
EOF

# --- STEP 3: AUTOMATED QC PLOTTING ---
echo "Merge complete. Launching R QC script for $TISSUE..."
Rscript /home/zw529/donglab/pipelines/scripts/QTL/_normQC.R "$TISSUE"

echo "Full QC pipeline complete for $TISSUE."
