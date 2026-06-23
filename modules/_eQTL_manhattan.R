#!/usr/bin/env Rscript
###########################################
# _eQTL_manhattan.R
# UPDATED: Generates both directional and gene-colored diagnostic plots
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

# Use FDR column if it exists, otherwise fall back to p-value
if ("FDR" %in% names(snp_plot)) {
    message("## Using FDR column for significance threshold")
    snp_plot[, log10p := -log10(`p-value`)]
    snp_plot[, fdr_val := FDR]
} else {
    message("## FDR column not found, using p-value")
    snp_plot[, log10p := -log10(`p-value`)]
    snp_plot[, fdr_val := `p-value`]
}

if ("FDR" %in% names(lead_plot)) {
    lead_plot[, log10p := -log10(`p-value`)]
    lead_plot[, fdr_val := FDR]
} else {
    lead_plot[, log10p := -log10(`p-value`)]
    lead_plot[, fdr_val := `p-value`]
}

sig_thresh_fdr <- fdr_cutoff

message("## Processing SNP coordinates (Autosomes 1-22 only)...")
prepare_coords <- function(df) {
    df[, chr_name := gsub("chr", "", as.character(chr), ignore.case = TRUE)]
    # Filter out Chromosome X or 23 explicitly
    df <- df[!(chr_name %in% c("X", "x", "23"))]
    df[, chr_num := as.integer(chr_name)]
    return(df[!is.na(chr_num)]) # Drops non-autosomes, Y, M, or parsing errors
}

snp_plot  <- prepare_coords(snp_plot)
lead_plot <- prepare_coords(lead_plot)

chr_info <- snp_plot[, .(chr_len = max(as.numeric(pos))), by = chr_num][order(chr_num)]
chr_info[, tot := shift(cumsum(as.numeric(chr_len)), fill = 0)]

setkey(snp_plot, chr_num); setkey(lead_plot, chr_num); setkey(chr_info, chr_num)

snp_plot  <- chr_info[snp_plot]
snp_plot[, cum_pos := pos + tot]

lead_plot <- chr_info[lead_plot]
lead_plot[, cum_pos := pos + tot]

snp_axis <- snp_plot[, .(centre = (max(cum_pos) + min(cum_pos)) / 2), by = chr_num]
extreme_hits <- lead_plot[order(log10p, decreasing = TRUE)][1:min(10, .N)]

# ========================================================================
# PLOT 1: ORIGINAL DIRECTIONAL COLORING (BETA-BASED)
# ========================================================================
message("## Assigning directional colors...")
snp_plot[, color_cat := as.character(chr_num %% 2)]
snp_plot[fdr_val <= sig_thresh_fdr & beta > 0, color_cat := "increase"]
snp_plot[fdr_val <= sig_thresh_fdr & beta < 0, color_cat := "decrease"]
snp_plot[telomeric_flag == TRUE, color_cat := "telomere"]

lead_plot[, color_cat := "none"]
lead_plot[fdr_val <= sig_thresh_fdr & beta > 0, color_cat := "inc_lead"]
lead_plot[fdr_val <= sig_thresh_fdr & beta < 0, color_cat := "dec_lead"]

message("## Creating Directional Plot...")
p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    geom_point(aes(colour=color_cat), size=0.6, alpha=0.5) +
    geom_point(data=lead_plot[color_cat != "none"], aes(x=cum_pos, y=log10p, colour=color_cat), shape=18, size=2.0) +
    geom_hline(yintercept=5, linetype="dashed", colour="grey30", linewidth=0.5) +
    geom_text_repel(data=extreme_hits, aes(x=cum_pos, y=log10p, label=geneid), size=2.5, colour="darkgreen", fontface="italic", box.padding = 0.5, max.overlaps = Inf) +
    scale_colour_manual(values=c("0"="#444444", "1"="#999999", "increase"="red", "decrease"="blue", "telomere"="#E69F00", "inc_lead"="red", "dec_lead"="blue", "none"="white")) +
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), breaks=snp_axis$centre, expand=c(0.01, 0.01)) +
    scale_y_continuous(expand=c(0.02, 0.5)) +
    labs(title=paste("cis-eQTL Manhattan Plot -", basename(out_prefix)), subtitle=paste("Reds: Increased expression | Blues: Decreased expression | FDR <", sig_thresh_fdr), x="Chromosome", y=expression(-log[10](p))) +
    theme_bw() + theme(legend.position="none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

out_file_png1 <- paste0(out_prefix, ".manhattan_by_SNP.png")
png(out_file_png1, width=18, height=6, units="in", res=300)
print(p1)
dev.off()

# ========================================================================
# PLOT 2: DIAGNOSTIC COLORING BY UNIQUE GENE ID (FIXED)
# ========================================================================
message("## Setting up Gene ID diagnostic colors...")

# Isolate significant versus background data frames to prevent scale conflicts
snp_sig  <- snp_plot[fdr_val <= sig_thresh_fdr]
snp_bg   <- snp_plot[fdr_val > sig_thresh_fdr]
lead_sig <- lead_plot[fdr_val <= sig_thresh_fdr]

message("## Creating Gene ID Diagnostic Plot...")
p2 := ggplot() +
    # 1. Background non-significant points (Zebra style using shape 21 allows independent 'fill')
    geom_point(data=snp_bg, aes(x=cum_pos, y=log10p, fill=as.character(chr_num %% 2)), 
               shape=21, stroke=0, size=0.6, alpha=0.3) +
    
    # 2. Significant points colored dynamically by their distinct Gene ID string
    geom_point(data=snp_sig, aes(x=cum_pos, y=log10p, colour=geneid), size=0.6, alpha=0.7) +
    
    # 3. Lead Diamonds mapped to the same Gene ID color scheme
    geom_point(data=lead_sig, aes(x=cum_pos, y=log10p, colour=geneid), shape=18, size=2.2) +
    
    geom_hline(yintercept=5, linetype="dashed", colour="grey30", linewidth=0.5) +
    geom_text_repel(data=extreme_hits, aes(x=cum_pos, y=log10p, label=geneid), 
                    size=2.5, colour="black", fontface="bold", box.padding = 0.5, max.overlaps = Inf) +
    
    # Control the background layers using fill, leaving scale_colour to handle genes dynamically
    scale_fill_manual(values=c("0"="#CCCCCC", "1"="#E5E5E5"), guide="none") +
    
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), breaks=snp_axis$centre, expand=c(0.01, 0.01)) +
    scale_y_continuous(expand=c(0.02, 0.5)) +
    
    labs(title=paste("Diagnostic eQTL Manhattan Plot (By Gene Locus) -", basename(out_prefix)),
         subtitle=sprintf("Every unique significant Gene ID is assigned an independent color | N Significant Genes = %d", uniqueN(snp_sig$geneid)),
         x="Chromosome", y=expression(-log[10](p))) +
    
    theme_bw() + 
    theme(legend.position="none", 
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())

out_file_png2 <- paste0(out_prefix, ".manhattan_by_GENE.png")
png(out_file_png2, width=18, height=6, units="in", res=300)
print(p2)
dev.off()

message(paste("## Directional plot saved to:", out_file_png1))
message(paste("## Gene diagnostic plot saved to:", out_file_png2))
message("## Done!")
