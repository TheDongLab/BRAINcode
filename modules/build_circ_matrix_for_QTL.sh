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

# Match paths
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
    circ_file = row["circ_path"]
    
    # DEFENSIVE CHECK: Skip if the path is a broken symlink or unreadable
    if not circ_file.is_file():
        print(f"Warning: Skipping {sample} - file missing or broken symlink: {circ_file}")
        continue
        
    try:
        circ = pd.read_csv(circ_file, sep=r"\s+")
    except Exception as e:
        print(f"Warning: Could not read {circ_file} due to error: {e}")
        continue
    
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

# Total Cohort Sums
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

# --- TRACK PER-SAMPLE PROFILE AVERAGES SAFELY ---
reads_filtered = reads_master[reads_master["circ_id"].isin(filtered_circ_ids)].copy()
mean_reads_per_circ = reads_filtered.drop(columns=["circ_id"]).mean(axis=1)

averages_summary = pd.DataFrame({
    "circ_id": reads_filtered["circ_id"],
    "avg_reads": mean_reads_per_circ
})

# Save Outputs
expr_final.to_csv(OUT / "circ_matrix.txt", sep="\t", index=False)
rna_found[["externalsampleid", "externalsubjectid", "tissue"]].to_csv(
    OUT / "circ_sample_metadata.csv", index=False
)

# Export individual data summaries for both plots
reads_summary[reads_summary["circ_id"].isin(filtered_circ_ids)].to_csv(OUT / "filtered_reads_summary.tmp", sep="\t", index=False)
averages_summary.to_csv(OUT / "filtered_averages_summary.tmp", sep="\t", index=False)

# GENERATE RUN STATISTICS
num_circs, num_samples = filtered_numeric.shape
total_cells = filtered_numeric.size
nonzero_cells = (filtered_numeric > 0).sum().sum()

print("\n" + "="*40)
print("        === MATRIX STATS ===")
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
print("        === CHROMOSOMES ===")
print("="*40)
chrs_series = filtered_circ_ids.reset_index(drop=True).str.split(':').str[0]
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
# Define global shared variables upfront to prevent inner-loop execution drops
y_ticks <- 10^(0:5)
out_dir <- "/home/zw529/donglab/data/target_ALS/QTL/"

# =========================================================================
# GRAPH 1: COHORT TOTAL DISTRIBUTION (SUMMED SPLIT-AXIS PLOT)
# =========================================================================
print("Generating original cohort total split-axis distribution plot and table...")

data_path <- paste0(out_dir, "filtered_reads_summary.tmp")
df_tot <- read.delim(data_path, header=TRUE, sep="\t")
df_tot <- df_tot[df_tot$total_reads > 0, ]

counts_table_tot <- table(df_tot$total_reads)
plot_data_tot <- data.frame(
    reads = as.numeric(names(counts_table_tot)),
    circ_count = as.numeric(counts_table_tot)
)

# --- SAVE TABLE 1: TOTALS DATA ---
write.table(plot_data_tot, paste0(out_dir, "circ_abundance_distribution_table.txt"), 
            sep="\t", row.names=FALSE, quote=FALSE)

df_tot_sorted <- df_tot[order(-df_tot$total_reads), ]
top_outliers_tot <- head(df_tot_sorted, 5)

png_out_tot <- paste0(out_dir, "circ_abundance_distribution.png")
png(png_out_tot, width=3400, height=1600, res=300)

layout(matrix(c(1, 2), nrow=1), widths=c(0.80, 0.20))

# Panel 1 (Totals 0 to 300)
par(mar=c(4.5, 6.5, 3, 0.5), bty="l") 
plot(plot_data_tot$reads, plot_data_tot$circ_count, log="y", type="h", col="red4", lwd=2, lend="square",
     xlim=c(0, 300), ylim=c(1, 100000), xlab="", ylab="", main="", yaxt="n")
axis(2, at=y_ticks, labels=format(y_ticks, big.mark=",", scientific=FALSE), las=1, cex.axis=1.0)
mtext("Number of circular RNAs", side=2, line=4.2, cex=1.1)
mtext("Number of back-spliced reads", side=1, line=2.5, at=185, cex=1.1)
mtext("Distribution of circRNA Expression by Back-spliced Read Support (Totals)", side=3, line=1, at=185, font=2, cex=1.2)

# Panel 2 (Totals 10k to 80k Outliers)
par(mar=c(4.5, 0.5, 3, 1.5), bty="n") 
plot(plot_data_tot$reads, plot_data_tot$circ_count, log="y", type="h", col="red4", lwd=2, lend="square",
     xlim=c(10000, 80000), ylim=c(1, 100000), xlab="", ylab="", main="", yaxt="n", xaxt="n")
axis(1, at=c(10000, 40000, 70000), labels=c("10k", "40k", "70k"))

par(xpd=TRUE)
for(i in 1:nrow(top_outliers_tot)) {
    x_pos <- top_outliers_tot$total_reads[i]
    if(x_pos >= 10000 & x_pos <= 80000) {
        segments(x0=x_pos, y0=10000, x1=x_pos, y1=3, lwd=1.2, col="black", lty=2)
        segments(x0=x_pos, y0=12000, x1=x_pos, y1=25000, lwd=1.5, col="black")
        text(x=x_pos, y=30000, labels=top_outliers_tot$circ_id[i], srt=90, adj=0, cex=0.5, font=2, col="black")
    }
}
text(2500, 0.5, "//", cex=1.2) 
dev.off()


# =========================================================================
# GRAPH 2: CROSS-SUBJECT AVERAGES DISTRIBUTION (AVERAGED SPLIT-AXIS PLOT)
# =========================================================================
print("Generating new cohort average split-axis distribution plot and table...")

avg_path <- paste0(out_dir, "filtered_averages_summary.tmp")
df_avg <- read.delim(avg_path, header=TRUE, sep="\t")
df_avg <- df_avg[df_avg$avg_reads > 0, ]

# Decreased bin width from 0.2 to 0.05 for sharper granularity
bin_width <- 0.05
breaks_seq <- seq(0, max(df_avg$avg_reads) + bin_width, by=bin_width)
h <- hist(df_avg$avg_reads, breaks=breaks_seq, plot=FALSE)

plot_data_avg <- data.frame(
    reads_bin_center = h$mids[h$counts > 0],
    circ_count = h$counts[h$counts > 0]
)

# --- SAVE TABLE 2: AVERAGES DATA ---
write.table(plot_data_avg, paste0(out_dir, "circ_abundance_averages_table.txt"), 
            sep="\t", row.names=FALSE, quote=FALSE)

df_avg_sorted <- df_avg[order(-df_avg$avg_reads), ]
top_outliers_avg <- head(df_avg_sorted, 5)

png_out_avg <- paste0(out_dir, "circ_abundance_averages.png")
png(png_out_avg, width=3400, height=1600, res=300)

layout(matrix(c(1, 2), nrow=1), widths=c(0.80, 0.20))

# Panel 1 (Averages 0 to 20) - Color updated to darkorange2
par(mar=c(4.5, 6.5, 3, 0.5), bty="l") 
plot(plot_data_avg$reads_bin_center, plot_data_avg$circ_count, log="y", type="h", col="darkorange2", lwd=1.5, lend="square",
     xlim=c(0, 20), ylim=c(1, 100000), xlab="", ylab="", main="", yaxt="n")
axis(2, at=y_ticks, labels=format(y_ticks, big.mark=",", scientific=FALSE), las=1, cex.axis=1.0)
mtext("Number of circular RNAs", side=2, line=4.2, cex=1.1)
mtext("Average number of back-spliced reads per sample", side=1, line=2.5, at=12.5, cex=1.1)
mtext("Distribution of circRNA Expression by Average Back-spliced Read Support", side=3, line=1, at=12.5, font=2, cex=1.2)

# Panel 2 (Averages Extreme Outliers) - Color updated to darkorange2
par(mar=c(4.5, 0.5, 3, 1.5), bty="n") 
max_avg <- max(plot_data_avg$reads_bin_center)
xlim_p2_avg <- c(max_avg * 0.15, max_avg * 1.05)

plot(plot_data_avg$reads_bin_center, plot_data_avg$circ_count, log="y", type="h", col="darkorange2", lwd=1.5, lend="square",
     xlim=xlim_p2_avg, ylim=c(1, 100000), xlab="", ylab="", main="", yaxt="n", xaxt="n")

axis(1, at=round(seq(xlim_p2_avg[1], xlim_p2_avg[2], length.out=3), 1))

par(xpd=TRUE)
for(i in 1:nrow(top_outliers_avg)) {
    x_pos <- top_outliers_avg$avg_reads[i]
    if(x_pos >= xlim_p2_avg[1] & x_pos <= xlim_p2_avg[2]) {
        segments(x0=x_pos, y0=10000, x1=x_pos, y1=3, lwd=1.2, col="black", lty=2)
        segments(x0=x_pos, y0=12000, x1=x_pos, y1=25000, lwd=1.5, col="black")
        # Text tilted to 45 degrees using srt=45, with adjusted base alignment
        text(x=x_pos, y=30000, labels=top_outliers_avg$circ_id[i], srt=45, adj=c(0, 0), cex=0.5, font=2, col="black")
    }
}

# Shifted the // break marker leftward out of the graph frame to prevent overlap
text(xlim_p2_avg[1] - (xlim_p2_avg[1] * 0.18), 0.5, "//", cex=1.2) 
dev.off()

print("SUCCESS: Both graphs (PNG) and both plotting coordinates (txt) saved cleanly.")
EOF
