#!/bin/bash
#SBATCH --job-name=build_and_plot_circ
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/build_and_plot_circ.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/build_and_plot_circ.err
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
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(scales)
})

out_dir <- "/home/zw529/donglab/data/target_ALS/QTL/"
table_path <- paste0(out_dir, "circ_abundance_distribution_table.txt")

print("Loading unified metrics table...")
df_metrics <- read.delim(table_path, header=TRUE, sep="\t")

y_ticks <- 10^(0:5)

# =========================================================================
# IDENTIFY TARGET GENES
# =========================================================================
als_genes <- c("ATXN1", "ATXN2", "HOMER1", "C9orf72", "SOD1", "FUS", "TARDBP", "TBK1")
syn_genes <- c("RIMS1", "RIMS2")
target_genes <- c(als_genes, syn_genes)

df_metrics$cat <- ifelse(grepl(paste(als_genes, collapse="|"), df_metrics$gene_name), "ALS",
                  ifelse(grepl(paste(syn_genes, collapse="|"), df_metrics$gene_name), "SYN", NA))

df_targets <- df_metrics[!is.na(df_metrics$cat) & df_metrics$total_reads >= 10, ]
cat_colors <- c("ALS" = "#E41A1C", "SYN" = "#377EB8")

print(paste("Found", nrow(df_targets), "target circRNAs with total reads >= 10 for annotation."))

# =========================================================================
# HELPER: save a ggplot/patchwork object as PNG + SVG (svglite backend -- no Cairo needed)
# =========================================================================
save_ggplot <- function(plot_obj, basename, width_in = 11.33, height_in = 5.33, dpi = 300) {
    ggsave(paste0(out_dir, basename, ".png"), plot=plot_obj, width=width_in, height=height_in, dpi=dpi, units="in")
    ggsave(paste0(out_dir, basename, ".svg"), plot=plot_obj, width=width_in, height=height_in, units="in", device=svglite::svglite)
}

# =========================================================================
# SHARED DATA: TOTAL-READS DISTRIBUTION
# =========================================================================
counts_table_tot <- table(df_metrics$total_reads)
plot_data_tot <- data.frame(
    reads = as.numeric(names(counts_table_tot)),
    circ_count = as.numeric(counts_table_tot)
)

df_tot_sorted <- df_metrics[order(-df_metrics$total_reads), ]
top_outliers_tot <- head(df_tot_sorted, 5)
outliers_p2_tot <- top_outliers_tot[top_outliers_tot$total_reads >= 10000 & top_outliers_tot$total_reads <= 80000, ]
if (nrow(outliers_p2_tot) > 0) {
    outliers_p2_tot$label_text <- paste0(outliers_p2_tot$circ_id, " (", outliers_p2_tot$gene_name, ")")
}

# Extra headroom past the rightmost outlier so rotated labels near the edge (e.g. circHOMER1)
# have room to render instead of getting clipped -- kept close to the original fixed limit,
# just nudged a touch higher rather than scaling dynamically off the data max.
xlim_p2_upper <- 85000

build_panel1_totals <- function(bar_color) {
    ggplot(plot_data_tot, aes(x=reads, y=circ_count)) +
        geom_segment(aes(xend=reads, y=0, yend=circ_count), color=bar_color, linewidth=1.1, lineend="square") +
        scale_y_log10(limits=c(1, 1e5), breaks=y_ticks, labels=comma) +
        coord_cartesian(xlim=c(0, 300)) +
        labs(x="Number of back-spliced reads", y="Number of circular RNAs",
             title="Distribution of circRNA Expression by Back-spliced Read Support (Totals)") +
        theme_classic(base_size=13) +
        theme(plot.title=element_text(size=12, face="bold"))
}

build_panel2_totals <- function(bar_color) {
    p <- ggplot(plot_data_tot, aes(x=reads, y=circ_count)) +
        geom_segment(aes(xend=reads, y=0, yend=circ_count), color=bar_color, linewidth=1.1, lineend="square") +
        scale_y_log10(limits=c(1, 1e5)) +
        scale_x_continuous(breaks=c(10000, 40000, 70000), labels=c("10k", "40k", "70k")) +
        coord_cartesian(xlim=c(10000, xlim_p2_upper), clip="off") +
        labs(x=NULL, y=NULL) +
        theme_classic(base_size=13) +
        theme(axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
              plot.margin=margin(t=5.5, r=35, b=5.5, l=5.5))

    if (nrow(outliers_p2_tot) > 0) {
        p <- p +
            geom_segment(data=outliers_p2_tot, aes(x=total_reads, xend=total_reads, y=3, yend=25000),
                         inherit.aes=FALSE, linewidth=0.5, linetype="dashed", color="black") +
            geom_text(data=outliers_p2_tot, aes(x=total_reads, y=30000, label=label_text),
                      inherit.aes=FALSE, angle=90, hjust=0, size=1.8, fontface="bold", color="black")
    }
    p
}

# =========================================================================
# GRAPH 1a: UNANNOTATED TOTAL DISTRIBUTION (plain red4 split-axis plot)
# =========================================================================
print("Generating unannotated cohort total split-axis distribution plot...")
plot_totals_unannotated <- build_panel1_totals("red4") + build_panel2_totals("red4") + plot_layout(widths=c(0.76, 0.24))
save_ggplot(plot_totals_unannotated, "circ_abundance_distribution_unannotated", width_in=12.5)

# =========================================================================
# GRAPH 1b: ANNOTATED TOTAL DISTRIBUTION (gray85, gene-labeled split-axis plot)
# =========================================================================
print("Generating annotated cohort total split-axis distribution plot...")

panel1_targets <- df_targets[df_targets$total_reads <= 300, ]
panel1_targets <- panel1_targets[order(panel1_targets$total_reads), ]
if (nrow(panel1_targets) > 0) {
    panel1_targets$y_base <- as.numeric(counts_table_tot[as.character(panel1_targets$total_reads)])
    idx <- seq_len(nrow(panel1_targets))
    stagger_factor <- ((idx - 1) %% 4) + 1
    panel1_targets$y_label_pos <- panel1_targets$y_base * (1.8 + (stagger_factor * 1.5))
    panel1_targets$rotation_angle <- c(45, 60, 75)[((idx - 1) %% 3) + 1]
    panel1_targets$y_arrow_end <- panel1_targets$y_base * 1.8
    panel1_targets$y_arrow_start <- panel1_targets$y_label_pos * 0.9
    panel1_targets$label_text <- paste0(
        "circ", regmatches(panel1_targets$gene_name, regexpr(paste(target_genes, collapse="|"), panel1_targets$gene_name)),
        " (", panel1_targets$total_reads, ")"
    )
}

panel2_targets <- df_targets[df_targets$total_reads >= 10000 & df_targets$total_reads <= 80000, ]
panel2_targets <- panel2_targets[order(panel2_targets$total_reads), ]
if (nrow(panel2_targets) > 0) {
    panel2_targets$y_base <- as.numeric(counts_table_tot[as.character(panel2_targets$total_reads)])
    panel2_targets$y_base[is.na(panel2_targets$y_base)] <- 1
    idx <- seq_len(nrow(panel2_targets))
    stagger_factor <- ((idx - 1) %% 4) + 1
    panel2_targets$y_label_pos <- 400 * (1.5 ^ stagger_factor)
    panel2_targets$rotation_angle <- c(35, 50, 65)[((idx - 1) %% 3) + 1]
    panel2_targets$y_arrow_end <- panel2_targets$y_base * 2.5
    panel2_targets$y_arrow_start <- panel2_targets$y_label_pos * 0.85
    panel2_targets$label_text <- paste0(
        "circ", regmatches(panel2_targets$gene_name, regexpr(paste(target_genes, collapse="|"), panel2_targets$gene_name)),
        " (", format(panel2_targets$total_reads, big.mark=","), ")"
    )
}

# Gene category colors are explained via a subtitle instead of a side legend
# (a legend was pushing in from the right and clipping the outlier panel).
cat_labels <- c("ALS" = "ALS-associated", "SYN" = "Synapse-associated")
fill_scale  <- scale_fill_manual(values=cat_colors, guide="none")
color_scale <- scale_color_manual(values=cat_colors, guide="none")

p1_annot <- build_panel1_totals("gray85") +
    labs(subtitle="Red = ALS-associated genes   |   Blue = Synapse-associated genes")
if (nrow(panel1_targets) > 0) {
    p1_annot <- p1_annot +
        geom_segment(data=panel1_targets, aes(x=total_reads, xend=total_reads, y=y_arrow_start, yend=y_arrow_end, color=cat),
                     inherit.aes=FALSE, linewidth=0.4, arrow=arrow(length=unit(0.06, "cm"), angle=15)) +
        geom_point(data=panel1_targets, aes(x=total_reads, y=y_base, fill=cat),
                   inherit.aes=FALSE, shape=21, color="black", size=2, stroke=0.3) +
        geom_text(data=panel1_targets, aes(x=total_reads, y=y_label_pos, label=label_text, angle=rotation_angle, color=cat),
                  inherit.aes=FALSE, hjust=0, size=1.9, fontface="italic") +
        fill_scale + color_scale
}

p2_annot <- build_panel2_totals("gray85")
if (nrow(panel2_targets) > 0) {
    p2_annot <- p2_annot +
        geom_segment(data=panel2_targets, aes(x=total_reads, xend=total_reads, y=y_arrow_start, yend=y_arrow_end, color=cat),
                     inherit.aes=FALSE, linewidth=0.4, arrow=arrow(length=unit(0.06, "cm"), angle=15)) +
        geom_point(data=panel2_targets, aes(x=total_reads, y=y_base, fill=cat),
                   inherit.aes=FALSE, shape=21, color="black", size=2, stroke=0.3) +
        geom_text(data=panel2_targets, aes(x=total_reads, y=y_label_pos, label=label_text, angle=rotation_angle, color=cat),
                  inherit.aes=FALSE, hjust=0, size=1.9, fontface="italic") +
        fill_scale + color_scale
}

plot_totals_annotated <- p1_annot + p2_annot + plot_layout(widths=c(0.76, 0.24))
save_ggplot(plot_totals_annotated, "circ_abundance_distribution", width_in=12.5)

# =========================================================================
# GRAPH 2: CROSS-SUBJECT AVERAGES DISTRIBUTION (unannotated, darkorange2)
# =========================================================================
print("Generating cohort average split-axis distribution plot...")

bin_width <- 0.05
breaks_seq <- seq(0, max(df_metrics$avg_reads) + bin_width, by=bin_width)
h <- hist(df_metrics$avg_reads, breaks=breaks_seq, plot=FALSE)

plot_data_avg <- data.frame(
    reads_bin_center = h$mids[h$counts > 0],
    circ_count = h$counts[h$counts > 0]
)

df_avg_sorted <- df_metrics[order(-df_metrics$avg_reads), ]
top_outliers_avg <- head(df_avg_sorted, 5)

max_avg <- max(plot_data_avg$reads_bin_center)
xlim_p2_avg <- c(max_avg * 0.15, max_avg * 1.05)

outliers_p2_avg <- top_outliers_avg[top_outliers_avg$avg_reads >= xlim_p2_avg[1] & top_outliers_avg$avg_reads <= xlim_p2_avg[2], ]
if (nrow(outliers_p2_avg) > 0) {
    outliers_p2_avg$label_text <- paste0(outliers_p2_avg$circ_id, " (", outliers_p2_avg$gene_name, ")")
}

p1_avg <- ggplot(plot_data_avg, aes(x=reads_bin_center, y=circ_count)) +
    geom_segment(aes(xend=reads_bin_center, y=0, yend=circ_count), color="darkorange2", linewidth=1, lineend="square") +
    scale_y_log10(limits=c(1, 1e5), breaks=y_ticks, labels=comma) +
    coord_cartesian(xlim=c(0, 20)) +
    labs(x="Average number of back-spliced reads per sample", y="Number of circular RNAs",
         title="Distribution of circRNA Expression by Average Back-spliced Read Support") +
    theme_classic(base_size=13) +
    theme(plot.title=element_text(size=12, face="bold"))

p2_avg <- ggplot(plot_data_avg, aes(x=reads_bin_center, y=circ_count)) +
    geom_segment(aes(xend=reads_bin_center, y=0, yend=circ_count), color="darkorange2", linewidth=1, lineend="square") +
    scale_y_log10(limits=c(1, 1e5)) +
    scale_x_continuous(breaks=round(seq(xlim_p2_avg[1], xlim_p2_avg[2], length.out=3), 1)) +
    coord_cartesian(xlim=xlim_p2_avg) +
    labs(x=NULL, y=NULL) +
    theme_classic(base_size=13) +
    theme(axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

if (nrow(outliers_p2_avg) > 0) {
    p2_avg <- p2_avg +
        geom_segment(data=outliers_p2_avg, aes(x=avg_reads, xend=avg_reads, y=3, yend=25000),
                     inherit.aes=FALSE, linewidth=0.5, linetype="dashed", color="black") +
        geom_text(data=outliers_p2_avg, aes(x=avg_reads, y=30000, label=label_text),
                  inherit.aes=FALSE, angle=45, hjust=0, size=1.8, fontface="bold", color="black")
}

plot_averages <- p1_avg + p2_avg + plot_layout(widths=c(0.80, 0.20))
save_ggplot(plot_averages, "circ_abundance_averages")

print("SUCCESS: circ matrix rebuilt; 3 plots (unannotated totals, annotated totals, averages) saved as PNG and SVG via ggplot2/ggsave.")
EOF
