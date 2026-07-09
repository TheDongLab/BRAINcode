#!/usr/bin/env Rscript
###########################################
# _cQTL_manhattan.R
# UPDATED: Generates both directional and circ-colored diagnostic plots
# ADDED: Optional interaction-mode baseline direction mapping to fix the red-wall bias
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

# Standardize column naming if necessary (MatrixEQTL uses gene, postprocess outputs circ_id)
if("gene" %in% names(snp_plot)) setnames(snp_plot, "gene", "circ_id")
if("gene" %in% names(lead_plot)) setnames(lead_plot, "gene", "circ_id")

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

message("## Processing SNP coordinates dynamically...")
prepare_coords <- function(df) {
    c_names <- gsub("chr", "", as.character(df$chr), ignore.case = TRUE)
    
    if (any(c_names %in% c("X", "x", "23")) && all(c_names %in% c("X", "x", "23", "NA", "", NA))) {
        df[, chr_num := 23]
    } else {
        df <- df[!(gsub("chr", "", as.character(chr), ignore.case = TRUE) %in% c("X", "x", "23"))]
        df[, chr_num := as.integer(gsub("chr", "", as.character(chr), ignore.case = TRUE))]
    }
    return(df[!is.na(chr_num)]) 
}

snp_plot  <- prepare_coords(snp_plot)
lead_plot <- prepare_coords(lead_plot)

# Determine if this is a dedicated standalone Chromosome X focused run
is_chrx_run <- all(snp_plot$chr_num == 23)

if (is_chrx_run) {
    message("## Setting up single-chromosome physical base-pair layout for ChrX...")
    snp_plot[, cum_pos := pos]
    lead_plot[, cum_pos := pos]
    
    max_pos <- max(snp_plot$pos, na.rm=TRUE)
    tick_breaks <- seq(0, max_pos, by = 25e6)
    tick_labels <- paste0(tick_breaks / 1e6, " Mb")
    x_expand <- c(0.02, 0.02)
    x_axis_label <- "Chromosome X Coordinate"
} else {
    message("## Setting up standard multi-chromosomal cumulative stacked layout...")
    chr_info <- snp_plot[, .(chr_len = max(as.numeric(pos))), by = chr_num][order(chr_num)]
    chr_info[, tot := shift(cumsum(as.numeric(chr_len)), fill = 0)]
    
    setkey(snp_plot, chr_num); setkey(lead_plot, chr_num); setkey(chr_info, chr_num)
    
    snp_plot  <- chr_info[snp_plot]
    snp_plot[, cum_pos := pos + tot]
    
    lead_plot <- chr_info[lead_plot]
    lead_plot[, cum_pos := pos + tot]
    
    snp_axis <- snp_plot[, .(centre = (max(cum_pos) + min(cum_pos)) / 2), by = chr_num]
    tick_breaks <- snp_axis$centre
    tick_labels <- paste0("chr", snp_axis$chr_num)
    x_expand <- c(0.01, 0.01)
    x_axis_label <- "Chromosome"
}

extreme_hits <- lead_plot[order(log10p, decreasing = TRUE)]
