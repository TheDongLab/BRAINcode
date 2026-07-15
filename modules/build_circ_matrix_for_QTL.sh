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
from collections import defaultdict

pd.set_option("future.no_silent_downcasting", True)

BASE = Path("/home/zw529/donglab/data/target_ALS")
OUT  = BASE / "QTL"
OUT.mkdir(exist_ok=True)

GENE_BED = Path("/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6")

def norm_id(x):
    return str(x).strip().replace("-", "_")

# =========================================================================
# STEP 1: LOAD GENOMIC ANNOTATIONS (BED6)
# =========================================================================
print("Loading GENCODE v49 gene annotations...")
genes_by_chrom = defaultdict(list)
if GENE_BED.exists():
    with open(GENE_BED, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                gene_name = parts[3]
                genes_by_chrom[chrom].append((start, end, gene_name))
    print(f"Loaded annotations for {len(genes_by_chrom)} chromosomes.")
else:
    print(f"Warning: Gene annotation BED6 file not found at {GENE_BED}. Table will use 'NA' for gene names.")

def find_overlapping_genes(circ_id):
    try:
        parts = circ_id.split(":")
        chrom = parts[0]
        coord_parts = parts[1].split("-")
        c_start = int(coord_parts[0])
        c_end = int(coord_parts[1])
    except Exception:
        return "NA"

    if chrom not in genes_by_chrom:
        return "NA"
    
    overlapping = []
    for g_start, g_end, g_name in genes_by_chrom[chrom]:
        if max(c_start, g_start) < min(c_end, g_end):
            overlapping.append(g_name)
            
    if overlapping:
        return ",".join(sorted(list(set(overlapping))))
    return "NA"

# =========================================================================
# STEP 2: MATRIX GENERATION & VARIANCE FILTERING
# =========================================================================
rna = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")
rna.columns = rna.columns.str.strip().str.lower()

all_circ = list(BASE.glob("**/RNAseq/Processed/*/circularRNA_known_circ_percentage.txt"))

def find_circ(sample_id):
    sid_norm = norm_id(sample_id)
    for c in all_circ:
        parent = norm_id(c.parent.name)
        if sid_norm == parent:
            return c
    return None

rna["circ_path"] = rna["externalsampleid"].apply(find_circ)
rna_found = rna[rna["circ_path"].notna()].copy()

rna_found = rna_found.drop_duplicates(subset=["externalsampleid"], keep="first")
print(f"Cleaned unique RNA samples for matrix processing: {len(rna_found)}")

expr = None
reads_list = []

for _, row in rna_found.iterrows():
    sample = norm_id(row["externalsampleid"])
    circ_file = row["circ_path"]
    
    if not circ_file.is_file():
        continue
        
    try:
        circ = pd.read_csv(circ_file, sep=r"\s+")
    except Exception:
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

circ_ids = expr['circ_id']
numeric_data = expr.drop('circ_id', axis=1)
row_sds = numeric_data.std(axis=1)

min_sd_threshold = 0.005
variant_mask = row_sds > min_sd_threshold

filtered_numeric = numeric_data[variant_mask].copy()
filtered_circ_ids = circ_ids[variant_mask]
expr_final = pd.concat([filtered_circ_ids.reset_index(drop=True), filtered_numeric.reset_index(drop=True)], axis=1)

# =========================================================================
# STEP 3: CALCULATE DETAILED ABUNDANCE METRICS (SINGLE UNIFIED TABLE)
# =========================================================================
print("Calculating read metrics and mapping host genes...")
reads_filtered = reads_master[reads_master["circ_id"].isin(filtered_circ_ids)].copy()
numeric_reads_filtered = reads_filtered.drop(columns=["circ_id"])

total_reads_filtered = numeric_reads_filtered.sum(axis=1)
mean_reads_filtered = numeric_reads_filtered.mean(axis=1)
mean_reads_nonzero_only = numeric_reads_filtered.replace(0, np.nan).mean(axis=1).fillna(0)

master_metrics = pd.DataFrame({
    "circ_id": reads_filtered["circ_id"],
    "total_reads": total_reads_filtered,
    "avg_reads": mean_reads_filtered,
    "avg_reads_nonzero_only": mean_reads_nonzero_only
})

master_metrics["gene_name"] = master_metrics["circ_id"].apply(find_overlapping_genes)

col_order = ["circ_id", "gene_name", "total_reads", "avg_reads", "avg_reads_nonzero_only"]
master_metrics = master_metrics[col_order]
master_metrics_filtered = master_metrics[master_metrics["total_reads"] > 0].copy()

expr_final.to_csv(OUT / "circ_matrix.txt", sep="\t", index=False)
rna_found[["externalsampleid", "externalsubjectid", "tissue"]].to_csv(
    OUT / "circ_sample_metadata.csv", index=False
)

master_metrics_filtered.to_csv(OUT / "circ_abundance_distribution_table.txt", sep="\t", index=False)

reads_summary = pd.DataFrame({"circ_id": reads_master["circ_id"], "total_reads": reads_master.drop(columns=["circ_id"]).sum(axis=1)})
reads_summary[reads_summary["circ_id"].isin(filtered_circ_ids)].to_csv(OUT / "filtered_reads_summary.tmp", sep="\t", index=False)
pd.DataFrame({"circ_id": reads_filtered["circ_id"], "avg_reads": mean_reads_filtered}).to_csv(OUT / "filtered_averages_summary.tmp", sep="\t", index=False)

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

# 2. Purge Python and load R cleanly to generate plots
module --force purge
module load R

Rscript - <<'EOF'
y_ticks <- 10^(0:5)
out_dir <- "/home/zw529/donglab/data/target_ALS/QTL/"
table_path <- paste0(out_dir, "circ_abundance_distribution_table.txt")

print("Loading unified metrics table...")
df_metrics <- read.delim(table_path, header=TRUE, sep="\t")

# =========================================================================
# IDENTIFY TARGET GENES FOR DYNAMIC PLOT ANNOTATION (>= 10 TOTAL READS)
# =========================================================================
target_genes <- c(
  "ATXN1", "ATXN2", "HOMER1", "QKI", "NTRK2", "SLC8A1", "RIMS1", "RIMS2", 
  "KCNN2", "MAP7", "FMN2", "SLAIN1", "MTCL1", "AKT3", "PTK2", "HIPK3", 
  "MIB1", "FBXW7", "VCAN", "LIFR", "RMST", "SOD1", "C9orf72", "FUS", 
  "TARDBP", "UNC13A", "KIF5A", "TBK1", "OPTN", "DLG4", "SYN1", "SNAP25", "SHANK3"
)

# Construct regex matching exact symbols surrounded by commas, triple underscores, or line boundaries
gene_regex <- paste0("(^|___|,)(", paste(target_genes, collapse="|"), ")($|___|,)")
df_targets <- df_metrics[grep(gene_regex, df_metrics$gene_name), ]
df_targets <- df_targets[df_targets$total_reads >= 10, ]

print(paste("Found", nrow(df_targets), "target circRNAs with total reads >= 10 for annotation."))

# =========================================================================
# GRAPH 1: COHORT TOTAL DISTRIBUTION (SUMMED SPLIT-AXIS PLOT)
# =========================================================================
print("Generating original cohort total split-axis distribution plot...")

counts_table_tot <- table(df_metrics$total_reads)
plot_data_tot <- data.frame(
    reads = as.numeric(names(counts_table_tot)),
    circ_count = as.numeric(counts_table_tot)
)

df_tot_sorted <- df_metrics[order(-df_metrics$total_reads), ]
top_outliers_tot <- head(df_tot_sorted, 5)

png_out_tot <- paste0(out_dir, "circ_abundance_distribution.png")
png(png_out_tot, width=3400, height=1600, res=300)

layout(matrix(c(1, 2), nrow=1), widths=c(0.80, 0.20))

# --- PANEL 1 (Totals 0 to 300) ---
par(mar=c(4.5, 6.5, 3, 0.5), bty="l") 
plot(plot_data_tot$reads, plot_data_tot$circ_count, log="y", type="h", col="red4", lwd=2, lend="square",
     xlim=c(0, 300), ylim=c(1, 100000), xlab="", ylab="", main="", yaxt="n")
axis(2, at=y_ticks, labels=format(y_ticks, big.mark=",", scientific=FALSE), las=1, cex.axis=1.0)
mtext("Number of circular RNAs", side=2, line=4.2, cex=1.1)
mtext("Number of back-spliced reads", side=1, line=2.5, at=185, cex=1.1)
mtext("Distribution of circRNA Expression by Back-spliced Read Support (Totals)", side=3, line=1, at=185, font=2, cex=1.2)

# Dynamically annotate target genes that fall in Panel 1 (10 to 300 reads)
par(xpd=TRUE)
panel1_targets <- df_targets[df_targets$total_reads >= 10 & df_targets$total_reads <= 300, ]
if(nrow(panel1_targets) > 0) {
    for(i in 1:nrow(panel1_targets)) {
        x_pos <- panel1_targets$total_reads[i]
        # Clean label to extract just the gene symbol name matching our list for space efficiency
        matched_gene <- regmatches(panel1_targets$gene_name[i], regexpr(paste(target_genes, collapse="|"), panel1_targets$gene_name[i]))
        label_text <- paste0(matched_gene, " (", x_pos, ")")
        
        # Plot an indicator point at the top of the bar, then draw a label above it
        points(x=x_pos, y=counts_table_tot[as.character(x_pos)], pch=18, col="blue", cex=0.8)
        text(x=x_pos, y=counts_table_tot[as.character(x_pos)] * 1.5, labels=label_text, srt=90, adj=0, cex=0.5, font=4, col="blue")
    }
}

# --- PANEL 2 (Totals 10k to 80k Outliers) ---
par(mar=c(4.5, 0.5, 3, 1.5), bty="n") 
plot(plot_data_tot$reads, plot_data_tot$circ_count, log="y", type="h", col="red4", lwd=2, lend="square",
     xlim=c(10000, 80000), ylim=c(1, 100000), xlab="", ylab="", main="", yaxt="n", xaxt="n")
axis(1, at=c(10000, 40000, 70000), labels=c("10k", "40k", "70k"))

# Standard Outlier Labels for Panel 2
for(i in 1:nrow(top_outliers_tot)) {
    x_pos <- top_outliers_tot$total_reads[i]
    label_text <- paste0(top_outliers_tot$circ_id[i], " (", top_outliers_tot$gene_name[i], ")")
    if(x_pos >= 10000 & x_pos <= 80000) {
        segments(x0=x_pos, y0=10000, x1=x_pos, y1=3, lwd=1.2, col="black", lty=2)
        segments(x0=x_pos, y0=12000, x1=x_pos, y1=25000, lwd=1.5, col="black")
        text(x=x_pos, y=30000, labels=label_text, srt=90, adj=0, cex=0.5, font=2, col="black")
    }
}

# Also check if any target genes were large enough to fall into Panel 2
panel2_targets <- df_targets[df_targets$total_reads >= 10000 & df_targets$total_reads <= 80000, ]
if(nrow(panel2_targets) > 0) {
    for(i in 1:nrow(panel2_targets)) {
        x_pos <- panel2_targets$total_reads[i]
        matched_gene <- regmatches(panel2_targets$gene_name[i], regexpr(paste(target_genes, collapse="|"), panel2_targets$gene_name[i]))
        text(x=x_pos, y=8000, labels=paste("★", matched_gene), col="blue", font=2, cex=0.6)
    }
}

text(2500, 0.5, "//", cex=1.2) 
dev.off()


# =========================================================================
# GRAPH 2: CROSS-SUBJECT AVERAGES DISTRIBUTION (AVERAGED SPLIT-AXIS PLOT)
# =========================================================================
print("Generating new cohort average split-axis distribution plot...")

bin_width <- 0.05
breaks_seq <- seq(0, max(df_metrics$avg_reads) + bin_width, by=bin_width)
h <- hist(df_metrics$avg_reads, breaks=breaks_seq, plot=FALSE)

plot_data_avg <- data.frame(
    reads_bin_center = h$mids[h$counts > 0],
    circ_count = h$counts[h$counts > 0]
)

df_avg_sorted <- df_metrics[order(-df_metrics$avg_reads), ]
top_outliers_avg <- head(df_avg_sorted, 5)

png_out_avg <- paste0(out_dir, "circ_abundance_averages.png")
png(png_out_avg, width=3400, height=1600, res=300)

layout(matrix(c(1, 2), nrow=1), widths=c(0.80, 0.20))

# Panel 1 (Averages 0 to 20)
par(mar=c(4.5, 6.5, 3, 0.5), bty="l") 
plot(plot_data_avg$reads_bin_center, plot_data_avg$circ_count, log="y", type="h", col="darkorange2", lwd=1.5, lend="square",
     xlim=c(0, 20), ylim=c(1, 100000), xlab="", ylab="", main="", yaxt="n")
axis(2, at=y_ticks, labels=format(y_ticks, big.mark=",", scientific=FALSE), las=1, cex.axis=1.0)
mtext("Number of circular RNAs", side=2, line=4.2, cex=1.1)
mtext("Average number of back-spliced reads per sample", side=1, line=2.5, at=12.5, cex=1.1)
mtext("Distribution of circRNA Expression by Average Back-spliced Read Support", side=3, line=1, at=12.5, font=2, cex=1.2)

# Panel 2 (Averages Extreme Outliers)
par(mar=c(4.5, 0.5, 3, 1.5), bty="n") 
max_avg <- max(plot_data_avg$reads_bin_center)
xlim_p2_avg <- c(max_avg * 0.15, max_avg * 1.05)

plot(plot_data_avg$reads_bin_center, plot_data_avg$circ_count, log="y", type="h", col="darkorange2", lwd=1.5, lend="square",
     xlim=xlim_p2_avg, ylim=c(1, 100000), xlab="", ylab="", main="", yaxt="n", xaxt="n")

axis(1, at=round(seq(xlim_p2_avg[1], xlim_p2_avg[2], length.out=3), 1))

par(xpd=TRUE)
for(i in 1:nrow(top_outliers_avg)) {
    x_pos <- top_outliers_avg$avg_reads[i]
    label_text <- paste0(top_outliers_avg$circ_id[i], " (", top_outliers_avg$gene_name[i], ")")
    if(x_pos >= xlim_p2_avg[1] & x_pos <= xlim_p2_avg[2]) {
        segments(x0=x_pos, y0=10000, x1=x_pos, y1=3, lwd=1.2, col="black", lty=2)
        segments(x0=x_pos, y0=12000, x1=x_pos, y1=25000, lwd=1.5, col="black")
        text(x=x_pos, y=30000, labels=label_text, srt=45, adj=c(0, 0), cex=0.5, font=2, col="black")
    }
}

text(xlim_p2_avg[1] - (xlim_p2_avg[1] * 0.18), 0.5, "//", cex=1.2) 
dev.off()

print("SUCCESS: Unified table saved, split-axis plots updated, and target genes annotated.")
EOF
