#!/usr/bin/env Rscript
###########################################
# _eQTL_manhattan.R
# Purpose: Generate two Manhattan plots from annotated cis-eQTL results:
#   1. By SNP position  (standard — where is the variant?)
#   2. By gene position (where is the affected gene?)
#
# Telomeric SNPs are plotted in a distinct colour but not removed.
# Lead SNPs (one per gene) are highlighted with a diamond shape.
#
# Usage:
#   Rscript _eQTL_manhattan.R \
#       annotated_eqtl.txt \
#       lead_snps.txt \
#       output_prefix \
#       [fdr_threshold=0.05]
###########################################

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
})

args       <- commandArgs(TRUE)
eqtl_file  <- args[1]   # full_annotated.txt
lead_file  <- args[2]   # lead_snps.txt
out_prefix <- args[3]
fdr_thresh <- ifelse(is.na(args[4]), 0.05, as.numeric(args[4]))

message("## Reading data...")
eqtl <- fread(eqtl_file)[chr %in% paste0("chr", 1:22)]
lead <- fread(lead_file)[chr %in% paste0("chr", 1:22)]

eqtl[, chr_num := as.integer(sub("chr","", chr))]
lead[, chr_num := as.integer(sub("chr","", chr))]

# Calculate p-value cutoff for the FDR line
fdr_pval_cutoff <- max(eqtl[FDR < fdr_thresh]$`p-value`, na.rm=TRUE)

# Function to build cumulative coordinates for the x-axis
build_cumulative <- function(dt, pos_col) {
    dt <- copy(dt)
    setorder(dt, chr_num, get(pos_col))
    chr_max <- dt[, .(max_pos = max(get(pos_col), na.rm=TRUE)), by=chr_num]
    chr_max[, offset := cumsum(shift(max_pos, fill=0))]
    dt <- merge(dt, chr_max[, .(chr_num, offset)], by="chr_num")
    dt[, cum_pos := get(pos_col) + offset]
    axis_df <- dt[, .(centre = mean(range(cum_pos, na.rm=TRUE))), by=chr_num]
    return(list(data=dt, axis=axis_df))
}

message("## Processing SNP coordinates...")
snp_build <- build_cumulative(eqtl[!is.na(pos)], "pos")
snp_data  <- snp_build$data
snp_axis  <- snp_build$axis

# Downsample non-significant points to keep PDF size manageable
set.seed(1234)
snp_sig    <- snp_data[`p-value` <= fdr_pval_cutoff]
snp_nonsig <- snp_data[`p-value` > fdr_pval_cutoff & `p-value` <= 0.001]
snp_nonsig <- snp_nonsig[sample(.N, min(.N, 150000))]
snp_plot   <- rbind(snp_sig, snp_nonsig)
snp_plot[, log10p := -log10(`p-value`)]

# Prepare lead SNPs for diamond shapes
lead_plot <- build_cumulative(lead[!is.na(pos)], "pos")$data
lead_plot[, log10p := -log10(`p-value`)]

# Identify extreme hits for "Needle" labels (-log10p > 100)
extreme_hits <- snp_plot[log10p > 100]

message("## Creating Plot...")
p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    # Normal and Telomeric points
    geom_point(data=snp_plot[!telomeric_flag], aes(colour=as.character(chr_num %% 2)), size=0.6, alpha=0.7) +
    geom_point(data=snp_plot[telomeric_flag == TRUE], colour="#E69F00", size=0.6, alpha=0.9) +
    # Lead SNPs (Diamonds)
    geom_point(data=lead_plot, aes(x=cum_pos, y=log10p), shape=18, size=2, colour="red") +
    # FDR Significance Line
    geom_hline(yintercept=-log10(fdr_pval_cutoff), linetype="dashed", colour="red", linewidth=0.5) +
    # NEEDLE LABELS for extreme hits
    geom_segment(data=extreme_hits, aes(x=cum_pos, xend=cum_pos, y=log10p, yend=log10p + 15), 
                 colour="darkgreen", linetype="dotted") +
    geom_text_repel(data=extreme_hits, aes(x=cum_pos, y=log10p + 15, label=geneid), 
                    size=3, colour="darkgreen", fontface="bold", box.padding = 0.5) +
    # Scales and Theme
    scale_colour_manual(values=c("0"="#555555", "1"="#AAAAAA")) +
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), breaks=snp_axis$centre, expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01, 0.1)) +
    labs(title="cis-eQTL Manhattan (SNP Position)", x="Chromosome", y=expression(-log[10](p))) +
    theme_bw() + 
    theme(legend.position="none", panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1))

out_file <- paste0(out_prefix, ".manhattan_by_SNP.pdf")
ggsave(out_file, p1, width=14, height=5)
message("## Done: ", out_file)
