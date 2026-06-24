#!/usr/bin/env python3
import gzip
import os
import sys

# Script location setup
GENCODE_DIR = os.path.dirname(os.path.abspath(__file__))
VCF_FILE = os.path.join(GENCODE_DIR, "GCF_000001405.40.gz")
OUTPUT_MAP = os.path.join(GENCODE_DIR, "ncbi_to_ucsc.txt")

def generate_map():
    if not os.path.exists(VCF_FILE):
        sys.exit(f"Error: Reference VCF file not found at {VCF_FILE}")
        
    print(f"Reading VCF header from: {os.path.basename(VCF_FILE)}")
    chrom_map = []
    autosome_count = 1
    
    with gzip.open(VCF_FILE, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                break  # Stop parsing immediately when header ends
                
            if line.startswith("##contig="):
                # Isolate the exact RefSeq accession ID
                parts = line.split("<ID=")[1].split(",")[0]
                ncbi_id = parts.strip()
                
                # Pair it sequentially with the standard UCSC naming scheme
                if autosome_count <= 22:
                    chrom_map.append((ncbi_id, f"chr{autosome_count}"))
                    autosome_count += 1
                elif autosome_count == 23:
                    chrom_map.append((ncbi_id, "chrX"))
                    autosome_count += 1
                elif autosome_count == 24:
                    chrom_map.append((ncbi_id, "chrY"))
                    autosome_count += 1
                elif autosome_count == 25:
                    chrom_map.append((ncbi_id, "chrM"))
                    autosome_count += 1

    # Write out the clean mapping file
    print(f"Writing mapping table to: {os.path.basename(OUTPUT_MAP)}")
    with open(OUTPUT_MAP, "w") as out:
        for ncbi, ucsc in chrom_map:
            out.write(f"{ncbi}\t{ucsc}\n")
            
    print("Mapping file successfully created.")

if __name__ == "__main__":
    generate_map()
