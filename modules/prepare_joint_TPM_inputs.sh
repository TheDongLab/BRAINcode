#!/bin/bash
#SBATCH --job-name=prep_joint_TPM
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/prep_joint_TPM.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/prep_joint_TPM.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

if [ $# -eq 0 ]; then
    echo "Usage: sbatch prepare_joint_TPM_inputs.sh <Tissue_Name>"
    echo "Example: sbatch prepare_joint_TPM_inputs.sh Motor_Cortex"
    exit 1
fi

TISSUE=$1
BASE="/home/zw529/donglab/data/target_ALS"
COV_FILE="$BASE/QTL/covariates.tsv"

# Path construction
TISSUE_DIR="$BASE/$TISSUE/RNAseq/Processed"
OUT_DIR="$BASE/$TISSUE/eQTL"
FINAL_TPM="$OUT_DIR/joint_TPM_matrix.tsv"

mkdir -p "$OUT_DIR"

module load Python/3.12.3-GCCcore-13.3.0

echo "Processing Tissue: $TISSUE"

# Find all normalization.tab files in the tissue-specific path
find "$TISSUE_DIR" -name "normalization.tab" ! -path "*/eQTL/*" > "$OUT_DIR/tpm_file_list.txt"

python3 - <<EOF
import pandas as pd
from pathlib import Path

# Config
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

# Load list of files
with open(list_file, 'r') as f:
    files = [line.strip() for line in f if line.strip()]

if not files:
    print(f"Error: No files found in $TISSUE_DIR")
    exit(1)

# Load metadata
meta = pd.read_csv(covariate_path, sep='\t')

# 1. Apply remapping to the metadata tissue column
meta['mapped_tissue'] = meta['tissue'].map(TISSUE_REMAP).fillna(meta['tissue'])

# 2. Filter for samples matching the target tissue (e.g. Motor_Cortex)
meta_subset = meta[meta['mapped_tissue'] == target_tissue].copy()

# 3. Handle ID normalization for matching
meta_subset['folder_id'] = meta_subset['externalsampleid'].str.replace('-', '_')
valid_ids = set(meta_subset['folder_id'])

joint_df = pd.DataFrame()
count = 0

for fpath in files:
    sample_id = Path(fpath).parent.name
    
    if sample_id in valid_ids:
        try:
            df = pd.read_csv(fpath, sep='\t', usecols=['gene_id', 'TPM'])
            df = df.rename(columns={'TPM': sample_id})
            
            if joint_df.empty:
                joint_df = df.set_index('gene_id')
            else:
                joint_df = joint_df.join(df.set_index('gene_id'), how='inner')
            count += 1
        except Exception as e:
            print(f"Error reading {sample_id}: {e}")

if not joint_df.empty:
    # Final name cleanup (underscores back to dashes)
    joint_df.columns = [c.replace('_', '-') for c in joint_df.columns]
    joint_df.to_csv(output_file, sep='\t')
    print(f"Successfully merged {count} samples into {output_file}")
    print(f"Matrix dimensions: {joint_df.shape}")
else:
    print("No samples matched the criteria.")
EOF
