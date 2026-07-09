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
    return(df[!is.na(chr_num)]) 
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

extreme_hits <- lead_plot[order(log10p, decreasing = TRUE)]

# ========================================================================
# GENERATE DIAGNOSTIC PLOTS (MANHATTAN BY GENE, BY DIRECTION & QQ)
# ========================================================================
out_pdf <- paste0(out_prefix, ".pdf")
pdf(out_pdf, width = 14, height = 5)

# 1. Setup Shared Metrics & Global Subsets
snp_plot$is_sig <- ifelse(snp_plot$fdr_val <= sig_thresh_fdr, "Significant", "Non-Significant")
snp_plot$bg_color <- ifelse(snp_plot$chr_num %% 2 == 0, "#7f8c8d", "#2c3e50")
label_hits <- head(extreme_hits[fdr_val <= sig_thresh_fdr], 10)
sig_thresh_line <- -log10(max(snp_plot$`p-value`[snp_plot$fdr_val <= sig_thresh_fdr], na.rm = TRUE))

# Shared theme parameters to minimize script bloating
manhattan_theme <- theme_bw() + theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none"
)

# ─────────────────────────────────────────────────────────────────
# LAYOUT 1: DIAGNOSTIC MANHATTAN PLOT (BY GENE LOCUS)
# ─────────────────────────────────────────────────────────────────
message("## Rendering Layout 1: Manhattan Overview by circRNA Locus...")

unique_sig_circs <- sort(unique(snp_plot$circ_id[snp_plot$is_sig == "Significant"]))
gene_palette <- scales::hue_pal()(length(unique_sig_circs))
names(gene_palette) <- unique_sig_circs

snp_plot$locus_color <- ifelse(snp_plot$is_sig == "Significant", snp_plot$circ_id, "Background")

p_locus <- ggplot(snp_plot, aes(x = cum_pos / 1e6, y = log10p)) +
    geom_point(data = subset(snp_plot, locus_color == "Background"), aes(colour = bg_color), size = 0.8, alpha = 0.4) +
    geom_point(data = subset(snp_plot, locus_color != "Background"), aes(colour = locus_color), size = 1.2, alpha = 0.8) +
    geom_hline(yintercept = sig_thresh_line, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text_repel(data = label_hits, aes(label = circ_id), size = 2.8, fontface = "bold", box.padding = 0.4, max.overlaps = 20) +
    scale_x_continuous(breaks = tick_breaks / 1e6, labels = tick_labels, expand = x_expand) +
    scale_colour_manual(values = c(gene_palette, "#7f8c8d" = "#D3D3D3", "#2c3e50" = "#A9A9A9")) +
    labs(title = paste("Diagnostic cQTL Manhattan Plot (By Gene Locus) -", basename(out_prefix)),
         subtitle = paste("N Significant Loci =", length(unique_sig_circs)), x = "Chromosome", y = expression(-log[10](p-value))) +
    manhattan_theme

print(p_locus)

# ─────────────────────────────────────────────────────────────────
# LAYOUT 2: DIRECTIONAL MANHATTAN PLOT (EFFECT DIRECTION)
# ─────────────────────────────────────────────────────────────────
message("## Rendering Layout 2: Manhattan Overview by Effect Direction...")

beta_col <- if("beta" %in% names(snp_plot)) "beta" else "slope"
snp_plot$direction <- ifelse(snp_plot$is_sig == "Non-Significant", "Background",
                             ifelse(snp_plot[[beta_col]] > 0, "Increased Expression", "Decreased Expression"))
snp_plot$dir_color <- ifelse(snp_plot$direction == "Background", snp_plot$bg_color, snp_plot$direction)

direction_palette <- c("Increased Expression" = "#E41A1C", "Decreased Expression" = "#377EB8", "#7f8c8d" = "#7f8c8d", "#2c3e50" = "#2c3e50")

p_direction <- ggplot(snp_plot, aes(x = cum_pos / 1e6, y = log10p)) +
    geom_point(data = subset(snp_plot, direction == "Background"), aes(colour = dir_color), size = 0.8, alpha = 0.5) +
    geom_point(data = subset(snp_plot, direction != "Background"), aes(colour = dir_color), size = 1.2, alpha = 0.8) +
    geom_hline(yintercept = sig_thresh_line, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text_repel(data = label_hits, aes(label = circ_id), size = 2.8, fontface = "bold", box.padding = 0.4, max.overlaps = 20, segment.color = "darkgreen") +
    scale_x_continuous(breaks = tick_breaks / 1e6, labels = tick_labels, expand = x_expand) +
    scale_colour_manual(values = direction_palette) +
    labs(title = paste("cis-cQTL Manhattan Plot -", basename(out_prefix)),
         subtitle = "Reds: Increased expression | Blues: Decreased expression", x = "Chromosome", y = expression(-log[10](p-value))) +
    manhattan_theme

print(p_direction)

# ─────────────────────────────────────────────────────────────────
# LAYOUT 3: THE Q-Q PLOT
# ─────────────────────────────────────────────────────────────────
message("## Rendering Layout 3: Q-Q Plot...")

observed_p <- sort(snp_plot$`p-value`)
qq_df <- data.frame(Expected = -log10((1:nrow(snp_plot)) / (nrow(snp_plot) + 1)), Observed = -log10(observed_p))

p_qq <- ggplot(qq_df, aes(x = Expected, y = Observed)) +
    geom_point(colour = "#34495e", size = 1.0, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, colour = "#e74c3c", linetype = "dashed", linewidth = 0.8) +
    labs(title = paste("Quantile-Quantile (Q-Q) Plot -", basename(out_prefix)),
         x = expression(Expected~~-log[10](p-value)), y = expression(Observed~~-log[10](p-value))) +
    theme_bw() + theme(panel.grid.minor = element_blank())

print(p_qq)
dev.off()

message("## Complete diagnostics pipeline complete! Output path: ", out_pdf)
