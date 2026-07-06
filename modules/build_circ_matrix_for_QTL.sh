#!/bin/bash
#SBATCH --job-name=build_circ_matrix
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/build_circ_matrix.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/build_circ_matrix.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=160G

set -euo pipefail

unset PYTHONPATH

# 1. Load Python for Matrix Generation & Filtering
module purge
module load Python/3.12.3-GCCcore-13.3.0
module load Python-bundle-PyPI/2024.06-GCCcore-13.3.0

python3 - <<'EOF'
import pandas as pd
import numpy as np
from pathlib import Path

pd.set_option("future.no_silent_downcasting", True)

BASE = Path("/home/zw529/donglab/data/target_ALS")
OUT  = BASE / "QTL"
OUT.mkdir(exist_ok=True)

def norm_id(x):
    return str(x).strip().replace("-", "_")

# Load Metadata
rna = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")
rna.columns = rna.columns.str.strip().str.lower()

# Find circularRNA files
all_circ = list(BASE.glob("**/RNAseq/Processed/*/circularRNA_known_circ_percentage.txt"))

def find_circ(sample_id):
    sid_norm = norm_id(sample_id)
    for c in all_circ:
        parent = norm_id(c.parent.name)
        if sid_norm == parent:
            return c
    return None

# Match paths (Garbage medical rows will naturally get None here since no folder matches them)
rna["circ_path"] = rna["externalsampleid"].apply(find_circ)
rna_found = rna[rna["circ_path"].notna()].copy()

# Drop true duplicated metadata entries to ensure clean 1:1 matrix headers
rna_found = rna_found.drop_duplicates(subset=["externalsampleid"], keep="first")
print(f"Cleaned unique RNA samples for matrix processing (CGND + non-CGND): {len(rna_found)}")

# Build Percentage Matrix & Track Raw Backspliced Reads
expr = None
reads_list = []

for _, row in rna_found.iterrows():
    sample = norm_id(row["externalsampleid"])
    circ = pd.read_csv(row["circ_path"], sep=r"\s+")
    
    if "circ_percent" not in circ.columns or "readNumber" not in circ.columns:
        continue

    circ["circ_id"] = (
        circ["chrom"].astype(str) + ":" + 
        circ["start"].astype(str) + "-" + 
        circ["end"].astype(str) + ":" + 
        circ["strand"].astype(str)
    )
    
    pct_df = circ[["circ_id", "circ_percent"]].copy()
    pct_df.columns = ["circ_id", sample]
    
    if expr is None:
        expr = pct_df
    else:
        expr = expr.merge(pct_df, on="circ_id", how="outer")
        
    read_df = circ[["circ_id", "readNumber"]].copy()
    read_df.columns = ["circ_id", sample]
    reads_list.append(read_df)

expr = expr.fillna(0).infer_objects(copy=False)

print("Aggregating cohort-wide back-spliced read counts...")
reads_master = reads_list[0]
for next_df in reads_list[1:]:
    reads_master = reads_master.merge(next_df, on="circ_id", how="outer")
reads_master = reads_master.fillna(0)

total_reads_per_circ = reads_master.drop(columns=["circ_id"]).sum(axis=1)
reads_summary = pd.DataFrame({"circ_id": reads_master["circ_id"], "total_reads": total_reads_per_circ})

# Variance-Based Filtering
circ_ids = expr['circ_id']
numeric_data = expr.drop('circ_id', axis=1)
row_sds = numeric_data.std(axis=1)

min_sd_threshold = 0.005
variant_mask = row_sds > min_sd_threshold

filtered_numeric = numeric_data[variant_mask].copy()
filtered_circ_ids = circ_ids[variant_mask]
expr_final = pd.concat([filtered_circ_ids.reset_index(drop=True), filtered_numeric.reset_index(drop=True)], axis=1)

# Save Matrix & Metadata Outputs
expr_final.to_csv(OUT / "circ_matrix.txt", sep="\t", index=False)
rna_found[["externalsampleid", "externalsubjectid", "tissue"]].to_csv(
    OUT / "circ_sample_metadata.csv", index=False
)

# Save intermediate reads summary for R plotting step
reads_summary[reads_summary["circ_id"].isin(filtered_circ_ids)].to_csv(OUT / "filtered_reads_summary.tmp", sep="\t", index=False)

# GENERATE RUN STATISTICS
num_circs, num_samples = filtered_numeric.shape
total_cells = filtered_numeric.size
nonzero_cells = (filtered_numeric > 0).sum().sum()

print("\n" + "="*40)
print("       === MATRIX STATS ===")
print("="*40)
print(f"Total circRNAs: {num_circs:,}")
print(f"Total Samples: {num_samples:,}")
print(f"Global Avg: {filtered_numeric.values.mean():.6f}")
if nonzero_cells > 0:
    print(f"Detected Avg (>0): {filtered_numeric.values[filtered_numeric.values > 0].mean():.6f}")
print(f"Sparsity: {(1 - (nonzero_cells / total_cells)) * 100:.2f}%")

active_counts = (filtered_numeric > 0).sum(axis=0)
top_sample = active_counts.idxmax()
print(f"Top Subject (>0%): {top_sample} ({active_counts[top_sample]:,} circs)")

print("\n" + "="*40)
print("       === CHROMOSOMES ===")
print("="*40)
chrs_series = filtered_circ_ids.str.split(':').str[0]
chr_counts = chrs_series.value_counts()
chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
for c in chr_order:
    if c in chr_counts:
        print(f"{c}: {chr_counts[c]:,}")
EOF

# 2. Purge Python and load R cleanly
module --force purge
module load R

Rscript - <<'EOF'
print("Generating abundance distribution plot using R...")

data_path <- "/home/zw529/donglab/data/target_ALS/QTL/filtered_reads_summary.tmp"
df <- read.delim(data_path, header=TRUE, sep="\t")

df <- df[df$total_reads > 0, ]

counts_table <- table(df$total_reads)
plot_data <- data.frame(
    reads = as.numeric(names(counts_table)),
    circ_count = as.numeric(counts_table)
)

png_out <- "/home/zw529/donglab/data/target_ALS/QTL/circ_abundance_distribution.png"
png(png_out, width=3000, height=1500, res=300)

par(mar=c(4.5, 4.5, 3, 1), bty="n")

plot(plot_data$reads, plot_data$circ_count, 
     log="y", 
     type="h", 
     col="red4", 
     lwd=2,
     lend="square",
     xlab="Number of back-spliced reads", 
     ylab="Number of circular RNAs",
     main="Distribution of circRNA Expression by Back-spliced Read Support",
     cex.lab=1.1, 
     cex.main=1.2,
     yaxt="n")

y_ticks <- 10^(0:5)
axis(2, at=y_ticks, labels=format(y_ticks, big.mark=","), las=1)

dev.off()

unlink(data_path)
print(paste("Plot saved successfully to:", png_out))
EOF
