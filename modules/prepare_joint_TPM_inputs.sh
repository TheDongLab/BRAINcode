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
VCF_FILE="$BASE/QTL/chromosome_joint_vcfs/target_ALS_chr22.vcf.gz"
OUT_DIR="$BASE/$TISSUE/eQTL"
FINAL_TPM="$OUT_DIR/joint_TPM_matrix.tsv"

mkdir -p "$OUT_DIR"

# --- Environment Setup ---
module purge
module load R/4.4.2-gfbf-2024a
module load Python/3.12.3-GCCcore-13.3.0
module load bcftools/1.20

echo "Processing Tissue: $TISSUE"

# --- STEP 1: GET GENOTYPE WHITELIST ---
echo "Extracting VCF sample IDs..."
VCF_IDS="$OUT_DIR/vcf_samples.txt"
bcftools query -l "$VCF_FILE" > "$VCF_IDS"

# --- STEP 2: DISCOVER ALL RNA-SEQ FILES ---
echo "Finding all RNA-seq files for $TISSUE..."
find "$BASE" -name "normalization.tab" | grep -i "$TISSUE" | grep "Processed" | grep -v "/eQTL/" > "$OUT_DIR/tpm_file_list.txt" || true

# --- STEP 3: PYTHON GENOTYPE-AWARE MERGE ---
python3 - <<EOF
import pandas as pd
from pathlib import Path
import sys

pd.set_option('future.no_silent_downcasting', True)

list_file = "$OUT_DIR/tpm_file_list.txt"
vcf_id_file = "$VCF_IDS"
covariate_path = "$COV_FILE"
output_file = "$FINAL_TPM"
target_tissue = "$TISSUE"

# 1. Load VCF Whitelist
with open(vcf_id_file, 'r') as f:
    vcf_set = set(line.strip() for line in f if line.strip())

# 2. Load Metadata (The Rosetta Stone)
meta = pd.read_csv(covariate_path, sep='\t')
meta.columns = meta.columns.str.strip()

# Note: Adjust column names if your metadata uses different headers for VCF vs RNA IDs
# Based on your previous data, 'externalsampleid' matches RNA folders. 
# We need the column that matches 'JHU1', 'GWF14-01', etc.
# Assuming 'subject_id' or similar contains the VCF-style name.
vcf_col = 'subject_id' 
rna_col = 'externalsampleid'

if vcf_col not in meta.columns:
    # Fallback: if subject_id doesn't exist, we look for a match in any column
    print(f"Warning: '{vcf_col}' not found. Attempting to find VCF ID column...")
    for col in meta.columns:
        if meta[col].isin(vcf_set).any():
            vcf_col = col
            print(f"Found VCF IDs in column: {vcf_col}")
            break

# 3. Create mapping: {RNA_ID: VCF_ID}
# Only keep entries where the VCF ID actually exists in our VCF file
mapping = meta[meta[vcf_col].isin(vcf_set)].set_index(rna_col)[vcf_col].to_dict()

# 4. Process files
try:
    with open(list_file, 'r') as f:
        files = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    sys.exit(1)

all_dfs = []
seen_vcf_ids = set()

for fpath in files:
    folder_name = Path(fpath).parent.name
    # Normalize folder name to match metadata (try underscore and dash)
    rna_id = None
    if folder_name in mapping:
        rna_id = folder_name
    elif folder_name.replace('_', '-') in mapping:
        rna_id = folder_name.replace('_', '-')
    
    if rna_id:
        vcf_id = mapping[rna_id]
        
        # Deduplication: Only take the first RNA-seq sample per Genotype ID
        if vcf_id in seen_vcf_ids:
            continue
            
        try:
            df = pd.read_csv(fpath, sep='\t', usecols=['gene_id', 'TPM'])
            df['gene_id'] = df['gene_id'].astype(str).str.strip().str.split('.').str[0]
            df = df.drop_duplicates(subset=['gene_id'])
            # CRITICAL: Rename column to VCF_ID so it matches Genotypes exactly
            df = df.set_index('gene_id').rename(columns={'TPM': vcf_id})
            all_dfs.append(df)
            seen_vcf_ids.add(vcf_id)
        except Exception as e:
            print(f"Error reading {folder_name}: {e}")

if all_dfs:
    print(f"Matched {len(all_dfs)} samples with genotypes. Merging...")
    joint_df = pd.concat(all_dfs, axis=1, join='outer').fillna(0).infer_objects(copy=False)
    
    # Filter: > 0.1 TPM in >= 20% of samples
    threshold = len(all_dfs) * 0.2
    mask = (joint_df > 0.1).sum(axis=1) >= threshold
    joint_df = joint_df[mask]
    
    joint_df.to_csv(output_file, sep='\t')
    print(f"Final Matrix dimensions: {joint_df.shape}")
    print(f"Used {len(seen_vcf_ids)} unique individuals.")
else:
    print("Zero matches between RNA-seq, Metadata, and VCF IDs.")
    sys.exit(1)
EOF

# --- STEP 4: AUTOMATED QC PLOTTING ---
echo "Merge complete. Launching R QC script for $TISSUE..."
Rscript /home/zw529/donglab/pipelines/scripts/QTL/_normQC.R "$TISSUE"

echo "Full QC pipeline complete for $TISSUE."
