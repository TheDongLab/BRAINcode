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

# --- THE FIX: STRICT LOGIC ORDER ---

# 1. Initialize EVERYTHING as Zebra (0 or 1) based on Chromosome
snp_plot[, color_cat := as.character(chr_num %% 2)]

# 2. ONLY overwrite the category IF significant AND non-telomeric
snp_plot[log10p >= sig_thresh & telomeric_flag == FALSE & beta > 0, color_cat := "increase"]
snp_plot[log10p >= sig_thresh & telomeric_flag == FALSE & beta < 0, color_cat := "decrease"]

# 3. Always apply telomere color last (overriding both zebra and significance)
snp_plot[telomeric_flag == TRUE, color_cat := "telomere"]

# 4. Lead SNPs (Diamonds) logic - strictly significant
lead_plot[, color_cat := "none"]
lead_plot[log10p >= sig_thresh & beta > 0, color_cat := "inc_lead"]
lead_plot[log10p >= sig_thresh & beta < 0, color_cat := "dec_lead"]

snp_axis <- snp_plot[, .(centre = (max(cum_pos) + min(cum_pos)) / 2), by = chr_num]
extreme_hits <- lead_plot[order(log10p, decreasing = TRUE)][1:min(10, .N)]

message("## Creating Plot...")

p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    # Background + Significant points
    geom_point(aes(colour=color_cat), size=0.7, alpha=0.6) +
    
    # Significant Lead SNPs (Diamonds)
    geom_point(data=lead_plot[color_cat != "none"], 
               aes(x=cum_pos, y=log10p, colour=color_cat), 
               shape=18, size=2.2) +
    
    geom_hline(yintercept=sig_thresh, linetype="dashed", colour="grey40", linewidth=0.5) +
    
    geom_text_repel(data=extreme_hits, 
                    aes(x=cum_pos, y=log10p, label=geneid), 
                    size=2.5, colour="black", fontface="italic", 
                    box.padding = 0.5, max.overlaps = Inf) +

    # Map the colors
    scale_colour_manual(values=c(
        "0"        = "#444444",   # Non-sig / Even Chr (Dark Grey)
        "1"        = "#999999",   # Non-sig / Odd Chr (Light Grey)
        "increase" = "red",       # Sig + Positive Beta
        "decrease" = "blue",      # Sig + Negative Beta
        "telomere" = "#E69F00",   # Orange
        "inc_lead" = "red",       # Lead Red Diamond
        "dec_lead" = "blue"       # Lead Blue Diamond
    )) +
    
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), 
                       breaks=snp_axis$centre, 
                       expand=c(0.01, 0.01)) +
    scale_y_continuous(expand=c(0.02, 0.5)) +
    
    labs(title=paste("cis-eQTL Manhattan Plot -", basename(out_prefix)),
         subtitle="Reds: Significant Increase | Blues: Significant Decrease | Greys: Non-significant",
         x="Chromosome", y=expression(-log[10](p))) +
    
    theme_bw() + 
    theme(legend.position="none", 
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())

out_file_png <- paste0(out_prefix, ".manhattan_by_SNP.png")
png(out_file_png, width=18, height=6, units="in", res=300)
print(p1)
dev.off()
