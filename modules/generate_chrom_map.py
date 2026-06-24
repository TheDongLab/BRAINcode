#!/usr/bin/env python3
import os

# Script location setup
GENCODE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_MAP = os.path.join(GENCODE_DIR, "ncbi_to_ucsc.txt")

def generate_map():
    # Official GRCh38.p14 / GCF_000001405.40 RefSeq to UCSC translation table
    # This covers chromosomes 1-22, X, Y, and Mitochondria exactly as defined by NCBI
    mapping_data = [
        ("NC_000001.11", "chr1"),
        ("NC_000002.12", "chr2"),
        ("NC_000003.12", "chr3"),
        ("NC_000004.12", "chr4"),
        ("NC_000005.10", "chr5"),
        ("NC_000006.12", "chr6"),
        ("NC_000007.14", "chr7"),
        ("NC_000008.11", "chr8"),
        ("NC_000009.12", "chr9"),
        ("NC_000010.11", "chr10"),
        ("NC_000011.10", "chr11"),
        ("NC_000012.12", "chr12"),
        ("NC_000013.11", "chr13"),
        ("NC_000014.9",  "chr14"),
        ("NC_000015.10", "chr15"),
        ("NC_000016.10", "chr16"),
        ("NC_000017.11", "chr17"),
        ("NC_000018.10", "chr18"),
        ("NC_000019.10", "chr19"),
        ("NC_000020.11", "chr20"),
        ("NC_000021.9",  "chr21"),
        ("NC_000022.11", "chr22"),
        ("NC_000023.11", "chrX"),
        ("NC_000024.10", "chrY"),
        ("NC_012920.1",  "chrM")  # GRCh38 mitochondrial genome accession
    ]
    
    print(f"Writing official GRCh38.p14 mapping table to: {os.path.basename(OUTPUT_MAP)}")
    
    with open(OUTPUT_MAP, "w") as out:
        for ncbi, ucsc in mapping_data:
            out.write(f"{ncbi}\t{ucsc}\n")
            
    print("Mapping file successfully created.")

if __name__ == "__main__":
    generate_map()
