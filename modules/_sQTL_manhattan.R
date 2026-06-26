#!/usr/bin/env Rscript
###########################################
# _sQTL_manhattan.R
###########################################
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(scales) # Required for advanced palette handling
})

args <- commandArgs(TRUE)
input_file  <- args[1]
lead_file   <- args[2]
out_prefix  <- args[3]
fdr_cutoff  <- as.numeric(args[4])

message("## Reading data...")
snp_plot  <- fread(input_file)
lead_plot <- fread(lead_file)

# Standardize feature column mapping for splicing
if("gene" %in% names(snp_plot)) setnames(snp_plot, "gene", "junction_id")
if("gene" %in% names(lead_plot)) setnames(lead_plot, "gene", "junction_id")

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

# Calculate cumulative positions
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

lead_plot[, color_cat := "none"]
lead_plot[fdr_val <= sig_thresh_fdr & beta > 0, color_cat := "inc_lead"]
lead_plot[fdr_val <= sig_thresh_fdr & beta < 0, color_cat := "dec_lead"]

message("## Creating Directional Plot...")
p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    geom_point(aes(colour=color_cat), size=0.6, alpha=0.5) +
    geom_point(data=lead_plot[color_cat != "none"], aes(x=cum_pos, y=log10p, colour=color_cat), shape=18, size=2.0) +
    geom_hline(yintercept=5, linetype="dashed", colour="grey30", linewidth=0.5) +
    geom_text_repel(data=extreme_hits, aes(x=cum_pos, y=log10p, label=junction_id), size=2.5, colour="darkgreen", fontface="italic", box.padding = 0.5, max.overlaps = Inf) +
    scale_colour_manual(values=c("0"="#444444", "1"="#999999", "increase"="red", "decrease"="blue", "inc_lead"="red", "dec_lead"="blue", "none"="white")) +
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), breaks=snp_axis$centre, expand=c(0.01, 0.01)) +
    scale_y_continuous(expand=c(0.02, 0.5)) +
    labs(title=paste("cis-sQTL Manhattan Plot -", basename(out_prefix)), subtitle=paste("Reds: Increased Intron Splicing (PSI) | Blues: Decreased Intron Splicing (PSI) | FDR <", sig_thresh_fdr), x="Chromosome", y=expression(-log[10](p))) +
    theme_bw() + theme(legend.position="none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

out_file_png1 <- paste0(out_prefix, ".manhattan_by_SNP.png")
png(out_file_png1, width=18, height=6, units="in", res=300)
print(p1)
dev.off()

# ========================================================================
# PLOT 2: HIGH-CONTRAST INTRON SEPARATION PALETTE (BY INTRON)
# ========================================================================
message("## Setting up high-contrast Intron ID diagnostic colors...")

snp_sig  <- snp_plot[fdr_val <= sig_thresh_fdr]
snp_bg   <- snp_plot[fdr_val > sig_thresh_fdr]
lead_sig <- lead_plot[fdr_val <= sig_thresh_fdr]

# Generate a high-contrast shuffled palette for significant introns
unique_introns <- unique(c(snp_sig$junction_id, lead_sig$junction_id))
n_introns <- length(unique_introns)

if (n_introns > 0) {
    # Generate default evenly spaced colors
    base_colors <- hue_pal()(n_introns)
    # Apply a deterministic prime-number skip shuffle to ensure neighboring introns get opposing color values
    shuffle_idx <- seq(1, n_introns, by = 1)
    if (n_introns > 5) {
        shuffle_idx <- order((shuffle_idx * 13) %% n_introns)
    }
    shuffled_colors <- base_colors[shuffle_idx]
    names(shuffled_colors) <- unique_introns
} else {
    shuffled_colors <- c()
}

message("## Creating Intron Diagnostic Plot...")
p2 <- ggplot() +
    # 1. Background non-significant points (Zebra style using shape 21)
    geom_point(data=snp_bg, aes(x=cum_pos, y=log10p, fill=as.character(chr_num %% 2)), 
               shape=21, stroke=0, size=0.6, alpha=0.3) +
    
    # 2. Significant points colored with shuffled contrasting palette
    geom_point(data=snp_sig, aes(x=cum_pos, y=log10p, colour=junction_id), size=0.6, alpha=0.8) +
    
    # 3. Lead Diamonds mapped to the same contrasting palette
    geom_point(data=lead_sig, aes(x=cum_pos, y=log10p, colour=junction_id), shape=18, size=2.2) +
    
    geom_hline(yintercept=5, linetype="dashed", colour="grey30", linewidth=0.5) +
    geom_text_repel(data=extreme_hits, aes(x=cum_pos, y=log10p, label=junction_id), 
                    size=2.5, colour="black", fontface="bold", box.padding = 0.5, max.overlaps = Inf) +
    
    # Background layers handled by fill, scale_colour handles the discrete contrasting colors
    scale_fill_manual(values=c("0"="#CCCCCC", "1"="#E5E5E5"), guide="none") +
    scale_colour_manual(values=shuffled_colors, guide="none") +
    
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), breaks=snp_axis$centre, expand=c(0.01, 0.01)) +
    scale_y_continuous(expand=c(0.02, 0.5)) +
    
    labs(title=paste("Diagnostic sQTL Manhattan Plot (By Splicing Intron) -", basename(out_prefix)),
         subtitle=sprintf("Neighboring overlapping introns are forced into highly distinct colors | N Significant Introns = %d", n_introns),
         x="Chromosome", y=expression(-log[10](p))) +
    
    theme_bw() + 
    theme(legend.position="none", 
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())

# REWORKED: Updated name destination mapping to intron target label
out_file_png2 <- paste0(out_prefix, ".manhattan_by_INTRON.png")
png(out_file_png2, width=18, height=6, units="in", res=300)
print(p2)
dev.off()

message(paste("## Directional plot saved to:", out_file_png1))
message(paste("## Intron diagnostic plot saved to:", out_file_png2))
message("## Done!")
