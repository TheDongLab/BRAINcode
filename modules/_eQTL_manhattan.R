#!/usr/bin/env Rscript
###########################################
# _eQTL_manhattan.R
# UPDATED: Directional coloring based on 'beta' using FDR column
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
    snp_plot[, log10p := -log10(`p-value`)]  # Still plot -log10(p) on y-axis
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

message(paste("## FDR threshold set to:", sig_thresh_fdr))
message(paste("## SNPs with FDR <", sig_thresh_fdr, "will be colored red/blue"))

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

message("## Assigning colors based on significance and beta direction...")

# --- CORRECTED COLOR ASSIGNMENT ---
# 1. Start EVERYONE as Zebra (Grey/Black based on chromosome)
snp_plot[, color_cat := as.character(chr_num %% 2)]

# 2. ONLY assign Red/Blue if BOTH conditions are met:
#    - FDR is significant AND
#    - beta has the correct direction
snp_plot[fdr_val <= sig_thresh_fdr & beta > 0, color_cat := "increase"]
snp_plot[fdr_val <= sig_thresh_fdr & beta < 0, color_cat := "decrease"]

# 3. Telomeres get priority (overrides everything)
snp_plot[telomeric_flag == TRUE, color_cat := "telomere"]

# 4. Lead SNPs (Diamonds) - same logic
lead_plot[, color_cat := "none"]
lead_plot[fdr_val <= sig_thresh_fdr & beta > 0, color_cat := "inc_lead"]
lead_plot[fdr_val <= sig_thresh_fdr & beta < 0, color_cat := "dec_lead"]

message("## Color Category Breakdown (SNPs):")
print(snp_plot[, .N, by=color_cat][order(color_cat)])

message("## Color Category Breakdown (Lead SNPs):")
print(lead_plot[, .N, by=color_cat][order(color_cat)])

snp_axis <- snp_plot[, .(centre = (max(cum_pos) + min(cum_pos)) / 2), by = chr_num]
extreme_hits <- lead_plot[order(log10p, decreasing = TRUE)][1:min(10, .N)]

message("## Creating Plot...")
p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    # Background + Significant points (Single layer)
    geom_point(aes(colour=color_cat), size=0.6, alpha=0.5) +
    
    # Significant Lead SNPs (Diamonds) - Plotted on TOP
    geom_point(data=lead_plot[color_cat != "none"], 
               aes(x=cum_pos, y=log10p, colour=color_cat), 
               shape=18, size=2.0) +
    
    # Threshold Line (visual reference for -log10(p) = 5)
    geom_hline(yintercept=5, linetype="dashed", colour="grey30", linewidth=0.5) +
    
    # Lead SNP labels
    geom_text_repel(data=extreme_hits, 
                    aes(x=cum_pos, y=log10p, label=geneid), 
                    size=2.5, colour="darkgreen", fontface="italic", 
                    box.padding = 0.5, max.overlaps = Inf) +
    
    # Color scheme
    scale_colour_manual(values=c(
        "0"        = "#444444",   # Even Chr Background (Dark)
        "1"        = "#999999",   # Odd Chr Background (Light)
        "increase" = "red",       # Sig ONLY (positive beta)
        "decrease" = "blue",      # Sig ONLY (negative beta)
        "telomere" = "#E69F00",   # Telomere Orange
        "inc_lead" = "red",       # Lead Diamond (positive beta)
        "dec_lead" = "blue",      # Lead Diamond (negative beta)
        "none"     = "white"      # Non-significant leads (shouldn't appear)
    )) +
    
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num), 
                       breaks=snp_axis$centre, expand=c(0.01, 0.01)) +
    scale_y_continuous(expand=c(0.02, 0.5)) +
    
    labs(title=paste("cis-eQTL Manhattan Plot -", basename(out_prefix)),
         subtitle=paste("Reds: Increased expression & FDR <", sig_thresh_fdr, "| Blues: Decreased expression & FDR <", sig_thresh_fdr,  "| Greys: Non-significant SNPs | Oranges: Telomeric SNPs"),
         x="Chromosome", y=expression(-log[10](p))) +
    
    theme_bw() + 
    theme(legend.position="none", 
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())

message("## Saving plot...")
out_file_png <- paste0(out_prefix, ".manhattan_by_SNP.png")
png(out_file_png, width=18, height=6, units="in", res=300)
print(p1)
dev.off()

message(paste("## Plot saved to:", out_file_png))
message("## Done!")
