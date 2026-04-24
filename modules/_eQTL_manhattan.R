#!/usr/bin/env Rscript
###########################################
# _eQTL_manhattan.R
# UPDATED: Directional coloring based on 'beta'
###########################################

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
})

args <- commandArgs(TRUE)
input_file  <- args[1]
lead_file   <- args[2]
out_prefix  <- args[3]
fdr_cutoff  <- as.numeric(args[4])

message("## Reading data...")
snp_plot  <- fread(input_file)
lead_plot <- fread(lead_file)

snp_plot[, log10p := -log10(`p-value`)]
lead_plot[, log10p := -log10(`p-value`)]

# Significance threshold for coloring
sig_thresh <- -log10(fdr_cutoff)

message("## Processing SNP coordinates...")
snp_plot[, chr_num := as.integer(gsub("chr", "", as.character(chr)))]
lead_plot[, chr_num := as.integer(gsub("chr", "", as.character(chr)))]

snp_plot <- snp_plot[!is.na(chr_num)]
lead_plot <- lead_plot[!is.na(chr_num)]

chr_info <- snp_plot[, .(chr_len = max(as.numeric(pos))), by = chr_num][order(chr_num)]
chr_info[, tot := shift(cumsum(as.numeric(chr_len)), fill = 0)]

setkey(snp_plot, chr_num); setkey(lead_plot, chr_num); setkey(chr_info, chr_num)

snp_plot <- chr_info[snp_plot]
snp_plot[, cum_pos := pos + tot]

lead_plot <- chr_info[lead_plot]
lead_plot[, cum_pos := pos + tot]

# Define direction only for Significant, Non-Telomeric SNPs
snp_plot[, color_cat := as.character(chr_num %% 2)] # Default to Zebra
snp_plot[log10p >= sig_thresh & beta > 0, color_cat := "increase"]
snp_plot[log10p >= sig_thresh & beta < 0, color_cat := "decrease"]
snp_plot[telomeric_flag == TRUE, color_cat := "telomere"]

# Apply same logic to Lead SNPs
lead_plot[, color_cat := "increase"]
lead_plot[beta < 0, color_cat := "decrease"]

snp_axis <- snp_plot[, .(centre = (max(cum_pos) + min(cum_pos)) / 2), by = chr_num]
extreme_hits <- lead_plot[order(log10p, decreasing = TRUE)][1:min(10, .N)]

message("## Creating Plot...")

p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    # 1. Plot all points using the combined color_cat
    # We use a single geom_point for the background/sig hits to keep it clean
    geom_point(aes(colour=color_cat), size=0.7, alpha=0.6) +
    
    # 2. Significant Lead SNPs get the Diamond shape
    geom_point(data=lead_plot[log10p >= sig_thresh], 
               aes(x=cum_pos, y=log10p, colour=color_cat), 
               shape=18, size=2.2) +
    
    geom_hline(yintercept=sig_thresh, linetype="dashed", colour="grey40", linewidth=0.5) +
    
    geom_text_repel(data=extreme_hits, 
                    aes(x=cum_pos, y=log10p, label=geneid), 
                    size=2.5, colour="black", fontface="italic", 
                    box.padding = 0.5, max.overlaps = Inf) +

    # 3. Scale Definition
    scale_colour_manual(values=c(
        "0" = "#444444",        # Even Chr Background
        "1" = "#999999",        # Odd Chr Background
        "telomere" = "#E69F00", # Orange
        "increase" = "#E41A1C", # Red
        "decrease" = "#377EB8"  # Blue
    )) +
    
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), 
                       breaks=snp_axis$centre, 
                       expand=c(0.01, 0.01)) +
    scale_y_continuous(expand=c(0.02, 0.5)) +
    
    labs(title=paste("cis-eQTL Manhattan Plot -", basename(out_prefix)),
         subtitle="Red: Increase Expr | Blue: Decrease Expr | Orange: Telomeric",
         x="Chromosome", y=expression(-log[10](p))) +
    
    theme_bw() + theme(legend.position="none", 
                       panel.grid.major.x = element_blank(),
                       panel.grid.minor = element_blank())

out_file_png <- paste0(out_prefix, ".manhattan_by_SNP.png")
png(out_file_png, width=18, height=6, units="in", res=300)
print(p1)
dev.off()
