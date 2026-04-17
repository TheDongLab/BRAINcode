#!/usr/bin/env python3
import pandas as pd
import sys
import os

def subset_matrix(matrix_path, samples_path, output_path):
    """
    Subsets a large expression matrix and reorders columns to match a sample list.
    """
    if not os.path.exists(samples_path):
        print(f"Error: Sample list {samples_path} not found.")
        sys.exit(1)

    # 1. Load the list of verified Matrix IDs (the 321 we found)
    with open(samples_path, 'r') as f:
        target_samples = [line.strip() for line in f if line.strip()]

    print(f"Targeting {len(target_samples)} samples for subsetting...")

    # 2. Check the matrix header to ensure all targets exist
    # Reading only the first row is extremely fast
    header = pd.read_csv(matrix_path, sep='\t', nrows=0)
    matrix_cols = header.columns.tolist()
    
    missing = [s for s in target_samples if s not in matrix_cols]
    if missing:
        print(f"Error: {len(missing)} samples missing from matrix: {missing[:5]}")
        sys.exit(1)

    # 3. Load only the necessary columns (Gene ID + our Targets)
    # This is the memory-efficient way to handle large genomic matrices
    print("Reading matrix columns...")
    try:
        # matrix_cols[0] is usually 'gene_id'
        df = pd.read_csv(matrix_path, sep='\t', usecols=[matrix_cols[0]] + target_samples, index_col=0)
    except Exception as e:
        print(f"Critical Error during load: {e}")
        sys.exit(1)

    # 4. Force the column order to match the sample list exactly
    # Matrix eQTL requires matrix columns to match covariate rows
    df = df[target_samples]

    # 5. Save the result
    df.to_csv(output_path, sep='\t')
    print(f"Successfully wrote subsetted matrix to: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python subset_matrix.py <input_matrix> <sample_list_file> <output_file>")
        sys.exit(1)
    subset_matrix(sys.argv[1], sys.argv[2], sys.argv[3])
