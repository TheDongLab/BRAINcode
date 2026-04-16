#!/usr/bin/env Rscript
###########################################
# _eQTL_manhattan.R
# FIXED: Sig-only red diamonds and chr labels
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

snp_axis <- snp_plot[, .(centre = (max(cum_pos) + min(cum_pos)) / 2), by = chr_num]
extreme_hits <- lead_plot[order(log10p, decreasing = TRUE)][1:min(10, .N)]

message("## Creating Plot...")

p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    # 1. Non-significant background SNPs (Zebra)
    geom_point(data=snp_plot[log10p < sig_thresh & telomeric_flag == FALSE], 
               aes(colour=as.character(chr_num %% 2)), 
               size=0.6, alpha=0.3) +
    
    # 2. Significant background SNPs (Still Zebra but brighter)
    geom_point(data=snp_plot[log10p >= sig_thresh & telomeric_flag == FALSE], 
               aes(colour=as.character(chr_num %% 2)), 
               size=0.7, alpha=0.6) +
    
    # 3. Telomeric SNPs
    geom_point(data=snp_plot[telomeric_flag == TRUE], 
               colour="#E69F00", size=0.6, alpha=0.8) +
    
    # 4. ONLY Significant Lead SNPs get the Red Diamond
    geom_point(data=lead_plot[log10p >= sig_thresh], 
               aes(x=cum_pos, y=log10p), 
               shape=18, size=1.8, colour="red") +
    
    geom_hline(yintercept=sig_thresh, linetype="dashed", colour="darkred", linewidth=0.5) +
    
    geom_text_repel(data=extreme_hits, 
                    aes(x=cum_pos, y=log10p, label=geneid), 
                    size=2.5, colour="darkgreen", fontface="italic", 
                    box.padding = 0.5, max.overlaps = Inf) +

    scale_colour_manual(values=c("0"="#444444", "1"="#999999")) +
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), 
                       breaks=snp_axis$centre, 
                       expand=c(0.01, 0.01)) +
    scale_y_continuous(expand=c(0.02, 0.5)) +
    
    labs(title=paste("cis-eQTL Manhattan Plot -", basename(out_prefix)),
         subtitle="Red Diamonds: Sig Lead SNPs | Orange: Telomeric SNPs",
         x="Chromosome", y=expression(-log[10](p))) +
    
    theme_bw() + theme(legend.position="none", panel.grid.major.x = element_blank())

out_file_png <- paste0(out_prefix, ".manhattan_by_SNP.png")
png(out_file_png, width=18, height=6, units="in", res=300)
print(p1)
dev.off()
