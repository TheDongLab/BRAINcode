#!/bin/bash
#SBATCH --job-name=prep_joint_TPM
#SBATCH --output=prep_joint_TPM_%j.out
#SBATCH --error=prep_joint_TPM_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

if [ $# -eq 0 ]; then
    echo "Usage: sbatch prepare_joint_TPM_inputs.sh <Tissue_Name>"
    echo "Example: sbatch prepare_joint_TPM_inputs.sh Cerebellum"
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
echo "Searching for files in: $TISSUE_DIR"

# Find all normalization.tab files in the tissue-specific path
find "$TISSUE_DIR" -name "normalization.tab" ! -path "*/eQTL/*" > "$OUT_DIR/tpm_file_list.txt"

python3 - <<EOF
import pandas as pd
from pathlib import Path
import sys

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
try:
    with open(list_file, 'r') as f:
        files = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print(f"Error: {list_file} not found.")
    sys.exit(1)

if not files:
    print(f"Error: No normalization.tab files found.")
    sys.exit(1)

# Load metadata
meta = pd.read_csv(covariate_path, sep='\t')

# 1. Apply remapping to the metadata tissue column
meta['mapped_tissue'] = meta['tissue'].map(TISSUE_REMAP).fillna(meta['tissue'])

# 2. Filter for samples matching the target tissue (e.g. Motor_Cortex)
meta_subset = meta[meta['mapped_tissue'] == target_tissue].copy()

# 3. Handle ID normalization for matching (dashes in meta -> underscores in folder names)
meta_subset['folder_id'] = meta_subset['externalsampleid'].str.replace('-', '_')
valid_ids = set(meta_subset['folder_id'])

all_dfs = []
print(f"Found {len(files)} total files. Filtering against {len(valid_ids)} valid metadata IDs...")

for fpath in files:
    # Folder name is the sample ID (e.g. CGND_HRA_00142)
    sample_id = Path(fpath).parent.name
    
    if sample_id in valid_ids:
        try:
            # Read only Gene ID and TPM
            df = pd.read_csv(fpath, sep='\t', usecols=['gene_id', 'TPM'])
            # Set index and rename TPM to the specific sample name
            df = df.set_index('gene_id').rename(columns={'TPM': sample_id})
            all_dfs.append(df)
        except Exception as e:
            print(f"Error reading {sample_id} at {fpath}: {e}")

if all_dfs:
    print(f"Merging {len(all_dfs)} dataframes...")
    # Concat along columns (axis=1)
    # join='inner' keeps only genes present in every single sample
    joint_df = pd.concat(all_dfs, axis=1, join='inner')
    
    # Final name cleanup (underscores back to dashes for R script compatibility)
    joint_df.columns = [c.replace('_', '-') for c in joint_df.columns]
    
    # Save the matrix
    joint_df.to_csv(output_file, sep='\t')
    print(f"Successfully merged {len(all_dfs)} samples into {output_file}")
    print(f"Final Matrix dimensions: {joint_df.shape}")
else:
    print("No samples matched the criteria. Check tissue names and folder IDs.")
EOF

echo "Process complete."
