#!/usr/bin/env Rscript
########################################################################
# _eQTL_regional_zoom.R
# Generates high-resolution coordinate-mapped single-chromosome zooms.
# Produces BOTH Directional (Beta-based) and Gene-Colored diagnostics.
########################################################################
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
})

args <- commandArgs(TRUE)
if (length(args) < 6) {
    stop("## ERROR: Missing arguments. Usage:\nRscript _eQTL_regional_zoom.R <cis_txt> <lead_snps_txt> <out_prefix> <fdr_cutoff> <target_chr> <start_mb> [end_mb]")
}

input_file  <- args[1]
lead_file   <- args[2]
out_prefix  <- args[3]
fdr_cutoff  <- as.numeric(args[4])
target_chr  <- gsub("chr", "", as.character(args[5]), ignore.case = TRUE)
start_mb    <- as.numeric(args[6])
end_mb      <- if(!is.na(args[7])) as.numeric(args[7]) else NULL

ref_snp_path <- "~/donglab/data/target_ALS/Cerebellum/eQTL/snp_location.txt"

# 1. Read Data safely
message("## Reading matrix tables...")
snp_raw  <- fread(input_file)
lead_raw <- fread(lead_file)

message("## Loading genomic position maps...")
if(!file.exists(ref_snp_path)) {
    stop(paste("## ERROR: Could not find reference location file at:", ref_snp_path))
}
snp_loc  <- fread(ref_snp_path, select = c("snpid", "chr", "pos"))

# Clean headers explicitly to prevent merge conflicts or suffix collision
setnames(snp_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(lead_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(snp_raw, "gene", "geneid", skip_absent = TRUE)
setnames(lead_raw, "gene", "geneid", skip_absent = TRUE)

snp_raw[, log10p := -log10(`p-value`)]
snp_raw[, fdr_val := if("FDR" %in% names(snp_raw)) FDR else `p-value`]

lead_raw[, log10p := -log10(`p-value`)]
lead_raw[, fdr_val := if("FDR" %in% names(lead_raw)) FDR else `p-value`]

# 2. Keyed Merge to avoid .x / .y suffix errors
message(paste0("## Filtering for Chromosome ", target_chr, "..."))
snp_loc[, chr_clean := gsub("chr", "", as.character(chr), ignore.case = TRUE)]
snp_loc_sub <- snp_loc[chr_clean == target_chr]

setkey(snp_raw, snpid)
setkey(lead_raw, snpid)
setkey(snp_loc_sub, snpid)

# Merge keeping coordinates explicitly clean
snp_zoom  <- snp_raw[snp_loc_sub, nomatch = NULL]
lead_zoom <- lead_raw[snp_loc_sub, nomatch = NULL]

# Window filters (converting Mb inputs to base-pairs)
if (!is.null(start_mb)) {
    snp_zoom  <- snp_zoom[pos >= (start_mb * 1e6)]
    lead_zoom <- lead_zoom[pos >= (start_mb * 1e6)]
}
if (!is.null(end_mb)) {
    snp_zoom  <- snp_zoom[pos <= (end_mb * 1e6)]
    lead_zoom <- lead_zoom[pos <= (end_mb * 1e6)]
}

if (nrow(snp_zoom) == 0) {
    stop("## ERROR: No data variants found matching that specific coordinate window.")
}

window_label <- if(!is.null(end_mb)) sprintf("chr%s:%.2f-%.2f Mb", target_chr, start_mb, end_mb) else sprintf("chr%s (>= %.2f Mb)", target_chr, start_mb)
extreme_hits <- lead_zoom[order(log10p, decreasing = TRUE)][1:min(6, .N)]

# ========================================================================
# PLOT A: DIRECTIONAL COLORING (BETA-BASED)
# ========================================================================
message("## Assigning directional colors for regional zoom...")
snp_zoom[, color_cat := "background"]
snp_zoom[fdr_val <= fdr_cutoff & beta > 0, color_cat := "increase"]
snp_zoom[fdr_val <= fdr_cutoff & beta < 0, color_cat := "decrease"]

lead_zoom[, color_cat := "none"]
lead_zoom[fdr_val <= fdr_cutoff & beta > 0, color_cat := "inc_lead"]
lead_zoom[fdr_val <= fdr_cutoff & beta < 0, color_cat := "dec_lead"]

p_dir <- ggplot() +
    geom_point(data=snp_zoom[color_cat == "background"], aes(x=pos / 1e6, y=log10p), colour="#999999", size=1.0, alpha=0.3) +
    geom_point(data=snp_zoom[color_cat == "increase"], aes(x=pos / 1e6, y=log10p), colour="red", size=1.4, alpha=0.7) +
    geom_point(data=snp_zoom[color_cat == "decrease"], aes(x=pos / 1e6, y=log10p), colour="blue", size=1.4, alpha=0.7) +
    geom_point(data=lead_zoom[color_cat == "inc_lead"], aes(x=pos / 1e6, y=log10p), colour="red", shape=18, size=3.5) +
    geom_point(data=lead_zoom[color_cat == "dec_lead"], aes(x=pos / 1e6, y=log10p), colour="blue", shape=18, size=3.5) +
    geom_hline(yintercept=5, linetype="dashed", colour="grey40", linewidth=0.5) +
    geom_text_repel(data=extreme_hits, aes(x=pos / 1e6, y=log10p, label=geneid), size=3.0, colour="black", fontface="bold.italic", box.padding = 0.6, max.overlaps = 15) +
    scale_x_continuous(expand=c(0.02, 0.02)) + scale_y_continuous(expand=c(0.02, 0.8)) +
    labs(title=paste("Directional Regional cis-eQTL Map -", basename(out_prefix)),
         subtitle=sprintf("Region: %s | Red: Increased Expression | Blue: Decreased Expression", window_label),
         x=paste("Chromosome", target_chr, "Position (Mb)"), y=expression(-log[10](p-value))) +
    theme_bw() + theme(panel.grid.minor = element_blank())

out_file_dir <- sprintf("%s.regional_zoom_dir_chr%s.png", out_prefix, target_chr)
png(out_file_dir, width=11, height=6, units="in", res=300)
print(p_dir)
dev.off()
message(paste("## Directional zoom plot saved to:", out_file_dir))

# ========================================================================
# PLOT B: DIAGNOSTIC COLORING BY UNIQUE GENE ID
# ========================================================================
message("## Setting up Gene ID diagnostic colors for regional zoom...")
snp_bg   <- snp_zoom[fdr_val > fdr_cutoff]
snp_sig  <- snp_zoom[fdr_val <= fdr_cutoff]
lead_sig <- lead_zoom[fdr_val <= fdr_cutoff]

p_gene <- ggplot() +
    geom_point(data=snp_bg, aes(x=pos / 1e6, y=log10p), colour="#999999", size=1.0, alpha=0.3) +
    geom_point(data=snp_sig, aes(x=pos / 1e6, y=log10p, colour=geneid), size=1.4, alpha=0.8) +
    geom_point(data=lead_sig, aes(x=pos / 1e6, y=log10p, colour=geneid), shape=18, size=3.5) +
    geom_hline(yintercept=5, linetype="dashed", colour="grey40", linewidth=0.5) +
    geom_text_repel(data=extreme_hits, aes(x=pos / 1e6, y=log10p, label=geneid), size=3.0, colour="black", fontface="bold", box.padding = 0.6, max.overlaps = 15) +
    scale_x_continuous(expand=c(0.02, 0.02)) + scale_y_continuous(expand=c(0.02, 0.8)) +
    labs(title=paste("Gene Diagnostic Regional cis-eQTL Map -", basename(out_prefix)),
         subtitle=sprintf("Region: %s | Every significant target gene structure gets a unique color", window_label),
         x=paste("Chromosome", target_chr, "Position (Mb)"), y=expression(-log[10](p-value))) +
    theme_bw() + theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size = 8.5), panel.grid.minor = element_blank())

out_file_gene <- sprintf("%s.regional_zoom_gene_chr%s.png", out_prefix, target_chr)
png(out_file_gene, width=11, height=6.5, units="in", res=300)
print(p_gene)
dev.off()
message(paste("## Gene diagnostic zoom plot saved to:", out_file_gene))
message("## All plots generated successfully.")
