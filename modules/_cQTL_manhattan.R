#!/usr/bin/env Rscript
###########################################
# _cQTL_manhattan.R
# UPDATED: Generates both directional and circ-colored diagnostic plots
# ADDED: Optional interaction-mode baseline direction mapping to fix the red-wall bias
# UPDATED: All plots strictly use FDR-based coloring for statistical significance
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

# Read interaction-specific optional arguments if provided
run_type   <- if(length(args) >= 5) args[5] else "standard"
std_path   <- if(length(args) >= 6) args[6] else NULL

message("## Reading data...")
snp_plot  <- fread(input_file)
lead_plot <- fread(lead_file)

# Standardize column naming for circular RNA structures
if("gene" %in% names(snp_plot)) setnames(snp_plot, "gene", "circ_id")
if("gene" %in% names(lead_plot)) setnames(lead_plot, "gene", "circ_id")

# Use FDR column if it exists, otherwise fall back to p-value
if ("FDR" %in% names(snp_plot)) {
    message("## Using FDR column for significance threshold")
    snp_plot$log10p <- -log10(snp_plot$`p-value`)
    snp_plot$fdr_val <- snp_plot$FDR
} else {
    message("## FDR column not found, using p-value")
    snp_plot$log10p <- -log10(snp_plot$`p-value`)
    snp_plot$fdr_val <- snp_plot$`p-value`
}

if ("FDR" %in% names(lead_plot)) {
    lead_plot$log10p <- -log10(lead_plot$`p-value`)
    lead_plot$fdr_val <- lead_plot$FDR
} else {
    lead_plot$log10p <- -log10(lead_plot$`p-value`)
    lead_plot$fdr_val <- lead_plot$`p-value`
}

sig_thresh_fdr <- fdr_cutoff

message("## Processing SNP coordinates dynamically...")
prepare_coords <- function(df) {
    c_names <- gsub("chr", "", as.character(df$chr), ignore.case = TRUE)
    
    if (any(c_names %in% c("X", "x", "23")) && all(c_names %in% c("X", "x", "23", "NA", "", NA))) {
        df$chr_num <- 23
    } else {
        df <- df[!(gsub("chr", "", as.character(chr), ignore.case = TRUE) %in% c("X", "x", "23"))]
        df$chr_num <- as.integer(gsub("chr", "", as.character(df$chr), ignore.case = TRUE))
    }
    return(df[!is.na(df$chr_num)]) 
}

snp_plot  <- prepare_coords(snp_plot)
lead_plot <- prepare_coords(lead_plot)

# Determine if this is a dedicated standalone Chromosome X focused run
is_chrx_run <- all(snp_plot$chr_num == 23)

if (is_chrx_run) {
    message("## Setting up single-chromosome physical base-pair layout for ChrX...")
    snp_plot$cum_pos <- snp_plot$pos
    lead_plot$cum_pos <- lead_plot$pos
    
    max_pos <- max(snp_plot$pos, na.rm=TRUE)
    tick_breaks <- seq(0, max_pos, by = 25e6)
    tick_labels <- paste0(tick_breaks / 1e6, " Mb")
    x_expand <- c(0.02, 0.02)
    x_axis_label <- "Chromosome X Coordinate"
} else {
    message("## Setting up standard multi-chromosomal cumulative stacked layout...")
    chr_info <- snp_plot[, .(chr_len = max(as.numeric(pos))), by = chr_num][order(chr_num)]
    chr_info$tot <- shift(cumsum(as.numeric(chr_info$chr_len)), fill = 0)
    
    setkey(snp_plot, chr_num); setkey(lead_plot, chr_num); setkey(chr_info, chr_num)
    
    snp_plot  <- chr_info[snp_plot]
    snp_plot$cum_pos <- snp_plot$pos + snp_plot$tot
    
    lead_plot <- chr_info[lead_plot]
    lead_plot$cum_pos <- lead_plot$pos + lead_plot$tot
    
    snp_axis <- snp_plot[, .(centre = (max(cum_pos) + min(cum_pos)) / 2), by = chr_num]
    tick_breaks <- snp_axis$centre
    tick_labels <- paste0("chr", snp_axis$chr_num)
    x_expand <- c(0.01, 0.01)
    x_axis_label <- "Chromosome"
}

extreme_hits <- lead_plot[order(log10p, decreasing = TRUE)][1:min(10, .N)]

# ========================================================================
# RE-MAP DIRECTIONALITY COLORS FOR INTERACTION MODES
# ========================================================================
eff_col <- if("slope" %in% names(snp_plot)) "slope" else "beta"

snp_plot$color_beta <- snp_plot[[eff_col]]
lead_plot$color_beta <- lead_plot[[eff_col]]

if (run_type == "interaction" && !is.null(std_path) && file.exists(std_path)) {
    message("## Interaction mode detected. Fetching baseline cQTL directionality for coloring...")
    
    std_df <- fread(std_path)
    if("gene" %in% names(std_df)) setnames(std_df, "gene", "circ_id")
    
    snp_plot$match_key <- paste0(snp_plot$snpid, "_", snp_plot$circ_id)
    lead_plot$match_key <- paste0(lead_plot$snpid, "_", lead_plot$circ_id)
    std_df$match_key <- paste0(std_df$snpid, "_", std_df$circ_id)
    
    std_betas <- setNames(std_df[[eff_col]], std_df$match_key)
    
    snp_plot$color_beta <- std_betas[snp_plot$match_key]
    lead_plot$color_beta <- std_betas[lead_plot$match_key]
    
    snp_plot$color_beta <- ifelse(is.na(snp_plot$color_beta), snp_plot[[eff_col]], snp_plot$color_beta)
    lead_plot$color_beta <- ifelse(is.na(lead_plot$color_beta), lead_plot[[eff_col]], lead_plot$color_beta)
}

# Set uniform FDR-based significance conditions and subtitle across modes
sig_cond_snp  <- snp_plot$fdr_val <= sig_thresh_fdr
sig_cond_lead <- lead_plot$fdr_val <= sig_thresh_fdr
sig_subtitle  <- paste("FDR <", sig_thresh_fdr)

# ========================================================================
# PLOT 1: DIRECTIONAL COLORING (MODIFIED SLOPE-BASED)
# ========================================================================
message("## Assigning directional colors...")
snp_plot$color_cat <- as.character(snp_plot$chr_num %% 2)

# Vectorized conditional mutations targeting directional changes
snp_plot$color_cat <- ifelse(sig_cond_snp & snp_plot$color_beta > 0, "increase", snp_plot$color_cat)
snp_plot$color_cat <- ifelse(sig_cond_snp & snp_plot$color_beta < 0, "decrease", snp_plot$color_cat)

if("telomeric_flag" %in% names(snp_plot)) {
    snp_plot$color_cat <- ifelse(snp_plot$telomeric_flag == TRUE, "telomere", snp_plot$color_cat)
}

lead_plot$color_cat <- "none"
lead_plot$color_cat <- ifelse(sig_cond_lead & lead_plot$color_beta > 0, "inc_lead", lead_plot$color_cat)
lead_plot$color_cat <- ifelse(sig_cond_lead & lead_plot$color_beta < 0, "dec_lead", lead_plot$color_cat)

message("## Creating Directional Plot...")
p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    geom_point(aes(colour=color_cat), size=0.6, alpha=0.5) +
    geom_point(data=subset(lead_plot, color_cat != "none"), aes(x=cum_pos, y=log10p, colour=color_cat), shape=18, size=2.0) +
    geom_hline(yintercept=5, linetype="dashed", colour="grey30", linewidth=0.5) +
    geom_text_repel(data=extreme_hits, aes(x=cum_pos, y=log10p, label=circ_id), size=2.5, colour="darkgreen", fontface="italic", box.padding = 0.5, max.overlaps = Inf) +
    scale_colour_manual(values=c("0"="#444444", "1"="#999999", "23"="#555555", "increase"="#E41A1C", "decrease"="#377EB8", "telomere"="#E69F00", "inc_lead"="#E41A1C", "dec_lead"="#377EB8", "none"="white")) +
    scale_x_continuous(labels=tick_labels, breaks=tick_breaks, expand=x_expand) +
    scale_y_continuous(breaks=c(0, 5, 10, 20, 30), expand=c(0.02, 0.5)) +
    labs(title=paste("cis-cQTL Manhattan Plot -", basename(out_prefix)), subtitle=paste("Reds: Increased expression | Blues: Decreased expression |", sig_subtitle), x=x_axis_label, y=expression(-log[10](p))) +
    theme_bw() + theme(legend.position="none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

out_file_png1 <- paste0(out_prefix, ".manhattan_by_SNP.png")
png(out_file_png1, width=18, height=6, units="in", res=300)
print(p1)
dev.off()

# ========================================================================
# PLOT 2: DIAGNOSTIC COLORING BY UNIQUE CIRC_ID (RANDOMIZED PALETTE)
# ========================================================================
message("## Setting up Circ ID diagnostic colors with a randomized palette...")

# Isolate significant versus background data frames strictly by FDR
snp_sig  <- snp_plot[fdr_val <= sig_thresh_fdr]
snp_bg   <- snp_plot[fdr_val > sig_thresh_fdr]
lead_sig <- lead_plot[fdr_val <= sig_thresh_fdr]

# 1. Get all unique significant IDs
unique_circs <- unique(snp_sig$circ_id)
n_circs      <- length(unique_circs)

if (n_circs > 0) {
    # 2. Generate a standard distinct palette matching the total count
    base_colors <- scales::hue_pal()(n_circs)
    
    # 3. Seed the randomizer so your plot colors remain reproducible across runs
    set.seed(100) 
    shuffled_colors <- sample(base_colors)
    
    # 4. Bind the shuffled colors explicitly to the unique IDs
    circ_palette <- setNames(shuffled_colors, unique_circs)
} else {
    circ_palette <- c()
}

message("## Creating Circ ID Diagnostic Plot...")
p2 <- ggplot() +
    geom_point(data=snp_bg, aes(x=cum_pos, y=log10p, fill=as.character(chr_num %% 2)), 
               shape=21, stroke=0, size=0.6, alpha=0.3) +
    geom_point(data=snp_sig, aes(x=cum_pos, y=log10p, colour=circ_id), size=0.6, alpha=0.7) +
    geom_point(data=lead_sig, aes(x=cum_pos, y=log10p, colour=circ_id), shape=18, size=2.2) +
    geom_hline(yintercept=5, linetype="dashed", colour="grey30", linewidth=0.5) +
    geom_text_repel(data=extreme_hits, aes(x=cum_pos, y=log10p, label=circ_id), 
                    size=2.5, colour="black", fontface="bold", box.padding = 0.5, max.overlaps = Inf) +
    scale_fill_manual(values=c("0"="#CCCCCC", "1"="#E5E5E5", "23"="#CCCCCC"), guide="none") +
    scale_colour_manual(values=circ_palette) +
    scale_x_continuous(labels=tick_labels, breaks=tick_breaks, expand=x_expand) +
    scale_y_continuous(breaks=c(0, 5, 10, 20, 30), expand=c(0.02, 0.5)) +
    labs(title=paste("Diagnostic cQTL Manhattan Plot (By Circ Locus) -", basename(out_prefix)),
         subtitle=sprintf("Unique significant Circ IDs are assigned a randomized color palette to maximize local contrast | N Significant Circs (%s) = %d", sig_subtitle, n_circs),
         x=x_axis_label, y=expression(-log[10](p))) +
    theme_bw() + 
    theme(legend.position="none", 
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())

out_file_png2 <- paste0(out_prefix, ".manhattan_by_CIRC.png")
png(out_file_png2, width=18, height=6, units="in", res=300)
print(p2)
dev.off()

message(paste("## Directional plot saved to:", out_file_png1))
message(paste("## Structural circ diagnostic plot saved to:", out_file_png2))
message("## Done!")
