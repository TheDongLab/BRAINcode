#!/usr/bin/env Rscript
###########################################
# _eQTL_manhattan.R
# Purpose: Generate a high-quality Manhattan plot from annotated eQTL results.
# Highlights: Lead SNPs, Telomeric SNPs, and top Gene Labels.
###########################################

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
})

args <- commandArgs(TRUE)
input_file  <- args[1] # .full_annotated.txt
lead_file   <- args[2] # .lead_snps.txt
out_prefix  <- args[3]
fdr_cutoff  <- as.numeric(args[4])

message("## Reading data...")
snp_plot  <- fread(input_file)
lead_plot <- fread(lead_file)

# Add log10 P-value column
snp_plot[, log10p := -log10(`p-value`)]
lead_plot[, log10p := -log10(`p-value`)]

message("## Processing SNP coordinates...")

# 1. Standardize Chromosome names to Integers for sorting
# Removes 'chr' prefix and converts to numeric (1-22)
snp_plot[, chr_num := as.integer(gsub("chr", "", as.character(chr)))]
lead_plot[, chr_num := as.integer(gsub("chr", "", as.character(chr)))]

# Remove any non-autosomal or unmapped SNPs if they exist
snp_plot <- snp_plot[!is.na(chr_num)]
lead_plot <- lead_plot[!is.na(chr_num)]

# 2. Calculate Cumulative Position (The Manhattan X-axis)
# We calculate the max position of each chromosome to stagger them
chr_info <- snp_plot[, .(chr_len = max(as.numeric(pos))), by = chr_num][order(chr_num)]
chr_info[, tot := shift(cumsum(as.numeric(chr_len)), fill = 0)]

# Merge the cumulative offsets back to the main data
setkey(snp_plot, chr_num)
setkey(lead_plot, chr_num)
setkey(chr_info, chr_num)

snp_plot <- chr_info[snp_plot]
snp_plot[, cum_pos := pos + tot]

lead_plot <- chr_info[lead_plot]
lead_plot[, cum_pos := pos + tot]

# 3. Create axis labels centered on each chromosome
snp_axis <- snp_plot[, .(centre = (max(cum_pos) + min(cum_pos)) / 2), by = chr_num]

# 4. Identify "Extreme Hits" for labeling (Top 10 most significant)
extreme_hits <- lead_plot[order(log10p, decreasing = TRUE)][1:min(10, .N)]

message("## Creating Plot...")

p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    # Background: Alternating colors for chromosomes (Zebra pattern)
    # Using explicit logical check for telomeric_flag
    geom_point(data=snp_plot[telomeric_flag == FALSE], 
               aes(colour=as.character(chr_num %% 2)), 
               size=0.6, alpha=0.5) +
    
    # Highlight Telomeric SNPs in Orange
    geom_point(data=snp_plot[telomeric_flag == TRUE], 
               colour="#E69F00", size=0.6, alpha=0.8) +
    
    # Highlight Lead SNPs as Red Diamonds
    geom_point(data=lead_plot, 
               aes(x=cum_pos, y=log10p), 
               shape=18, size=1.8, colour="red") +
    
    # Significance Threshold Line
    geom_hline(yintercept=-log10(fdr_cutoff), 
               linetype="dashed", colour="red", linewidth=0.5) +
    
    # Needles and Labels for the very top hits
    geom_segment(data=extreme_hits, 
                 aes(x=cum_pos, xend=cum_pos, y=0, yend=log10p), 
                 colour="darkgreen", linetype="dotted", linewidth=0.3) +
    
    geom_text_repel(data=extreme_hits, 
                    aes(x=cum_pos, y=log10p, label=geneid), 
                    size=2.5, colour="darkgreen", fontface="italic", 
                    box.padding = 0.5, max.overlaps = Inf) +

    # Color definitions (Dark/Light Gray)
    scale_colour_manual(values=c("0"="#444444", "1"="#999999")) +
    
    # Axis formatting
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), 
                       breaks=snp_axis$centre, 
                       expand=c(0.01, 0.01)) +
    scale_y_continuous(expand=c(0.02, 0.5)) +
    
    labs(title=paste("cis-eQTL Manhattan Plot -", basename(out_prefix)),
         subtitle="Red: Lead SNPs | Orange: Telomeric SNPs (<500kb)",
         x="Chromosome", 
         y=expression(-log[10](p))) +
    
    theme_bw() + 
    theme(
        legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size=8),
        plot.title = element_text(face="bold")
    )

# Using png() allows the cluster to "flatten" the 7 million dots into pixels
out_file_png <- paste0(out_prefix, ".manhattan_by_SNP.png")

# 300 DPI is high-enough quality for publication/presentation
png(out_file_png, width=18, height=6, units="in", res=300)
print(p1)
dev.off()

message("## Done: ", out_file_png)
