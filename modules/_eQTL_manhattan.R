#!/usr/bin/env Rscript
###########################################
# _eQTL_manhattan.R
# Purpose: Generate two Manhattan plots from annotated cis-eQTL results:
#   1. By SNP position  (standard — where is the variant?)
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
eqtl_file  <- args[1]
lead_file  <- args[2]
out_prefix <- args[3]
fdr_thresh <- ifelse(is.na(args[4]), 0.05, as.numeric(args[4]))

message("## Reading data...")
eqtl <- fread(eqtl_file)[chr %in% paste0("chr", 1:22)]
lead <- fread(lead_file)[chr %in% paste0("chr", 1:22)]

eqtl[, chr_num := as.integer(sub("chr","", chr))]
lead[, chr_num := as.integer(sub("chr","", chr))]

fdr_pval_cutoff <- max(eqtl[FDR < fdr_thresh]$`p-value`, na.rm=TRUE)

# Fixed build_cumulative function
build_cumulative <- function(dt, pos_col) {
    dt <- copy(dt)
    # Use setorderv for dynamic column names
    setorderv(dt, c("chr_num", pos_col))
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

set.seed(42)
snp_sig    <- snp_data[`p-value` <= fdr_pval_cutoff]
snp_nonsig <- snp_data[`p-value` > fdr_pval_cutoff & `p-value` <= 0.001]
snp_nonsig <- snp_nonsig[sample(.N, min(.N, 150000))]
snp_plot   <- rbind(snp_sig, snp_nonsig)
snp_plot[, log10p := -log10(`p-value`)]

lead_plot <- build_cumulative(lead[!is.na(pos)], "pos")$data
lead_plot[, log10p := -log10(`p-value`)]

# Extreme hits for "Needle" labels (-log10p > 100)
extreme_hits <- snp_plot[log10p > 100]

message("## Creating Plot...")
p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    # Use explicit logical checks to satisfy data.table scoping
    geom_point(data=snp_plot[telomeric_flag == FALSE], aes(colour=as.character(chr_num %% 2)), size=0.6, alpha=0.7) +
    geom_point(data=snp_plot[telomeric_flag == TRUE], colour="#E69F00", size=0.6, alpha=0.9) +
    
    # Lead SNPs (Diamonds)
    geom_point(data=lead_plot, aes(x=cum_pos, y=log10p), shape=18, size=2, colour="red") +
    geom_hline(yintercept=-log10(fdr_pval_cutoff), linetype="dashed", colour="red", linewidth=0.5) +
    
    # Needles starting from y=0
    geom_segment(data=extreme_hits, aes(x=cum_pos, xend=cum_pos, y=0, yend=log10p + 5), 
                 colour="darkgreen", linetype="dotted", linewidth=0.4) +
    geom_text_repel(data=extreme_hits, aes(x=cum_pos, y=log10p + 5, label=geneid), 
                    size=3, colour="darkgreen", fontface="italic", box.padding = 0.5) +

    # Scales and Theme
    scale_colour_manual(values=c("0"="#555555", "1"="#AAAAAA")) +
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), breaks=snp_axis$centre, expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01, 0.1)) +
    labs(title="cis-eQTL Manhattan (SNP Position)", x="Chromosome", y=expression(-log[10](p))) +
theme_bw() + 
    theme(
        legend.position="none", 
        panel.grid=element_blank(), 
        axis.text.x=element_text(angle=45, hjust=1, size=8) # Smaller text to fit all chrs
    )

out_file <- paste0(out_prefix, ".manhattan_by_SNP.pdf")
# WIDER CANVAS: 18x5 ensures all 22 chromosomes fit comfortably
ggsave(out_file, p1, width=18, height=5, limitsize=FALSE)
message("## Done: ", out_file)
