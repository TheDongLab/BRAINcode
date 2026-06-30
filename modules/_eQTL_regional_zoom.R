#!/usr/bin/env Rscript
########################################################################
# _eQTL_regional_zoom.R
# Dynamically zooms to coordinates OR anchors on a specific target SNP.
# Generates BOTH Directional (Beta-based) and Gene-Colored diagnostics.
########################################################################
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
})

args <- commandArgs(TRUE)
if (length(args) < 5) {
    stop("## ERROR: Missing arguments.\nUsage Option 1 (By Target SNP):  Rscript _eQTL_regional_zoom.R <cis_txt> <lead_snps_txt> <out_prefix> <fdr_cutoff> <target_snp>\nUsage Option 2 (By Region Below): Rscript _eQTL_regional_zoom.R <cis_txt> <lead_snps_txt> <out_prefix> <fdr_cutoff> <chr> <start_mb> <end_mb>")
}

input_file  <- args[1]
lead_file   <- args[2]
out_prefix  <- args[3]
fdr_cutoff  <- as.numeric(args[4])
query_param <- args[5] # Can be a SNP identifier (rs...) or a chromosome string

ref_snp_path <- "~/donglab/data/target_ALS/Cerebellum/eQTL/snp_location.txt"

# 1. Read Data
message("## Reading matrix tables...")
snp_raw  <- fread(input_file)
lead_raw <- fread(lead_file)

message("## Loading genomic position maps...")
if(!file.exists(ref_snp_path)) {
    stop(paste("## ERROR: Could not find reference location file at:", ref_snp_path))
}
snp_loc  <- fread(ref_snp_path, select = c("snpid", "chr", "pos"))

# Clean headers explicitly to align keys safely
setnames(snp_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(lead_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(snp_raw, "gene", "geneid", skip_absent = TRUE)
setnames(lead_raw, "gene", "geneid", skip_absent = TRUE)

snp_raw[, log10p := -log10(`p-value`)]
snp_raw[, fdr_val := if("FDR" %in% names(snp_raw)) FDR else `p-value`]
lead_raw[, log10p := -log10(`p-value`)]
lead_raw[, fdr_val := if("FDR" %in% names(lead_raw)) FDR else `p-value`]

# 2. Determine Coordinates dynamically based on parameter style
target_chr <- NULL; start_base <- NULL; end_base <- NULL; anchor_snp <- NULL

if (grepl("^rs", query_param, ignore.case = TRUE)) {
    anchor_snp <- query_param
    message(paste("## Look-up mode active: Finding anchors for SNP", anchor_snp))
    snp_lookup <- snp_loc[snpid == anchor_snp]
    if (nrow(snp_lookup) == 0) {
        stop(paste("## ERROR: Could not find coordinates for target SNP:", anchor_snp))
    }
    target_chr <- gsub("chr", "", as.character(snp_lookup$chr[1]), ignore.case = TRUE)
    snp_pos    <- as.numeric(snp_lookup$pos[1])
    # Build a tight 1 Mb locus window centered on the variant (+/- 500 kb)
    start_base <- max(0, snp_pos - 500000)
    end_base   <- snp_pos + 500000
    window_label <- sprintf("Locus window around %s (chr%s:%.2f-%.2f Mb)", anchor_snp, target_chr, start_base/1e6, end_base/1e6)
    out_suffix   <- paste0("_zoom_", anchor_snp)
} else {
    target_chr   <- gsub("chr", "", as.character(query_param), ignore.case = TRUE)
    start_base   <- as.numeric(args[6]) * 1e6
    end_base     <- as.numeric(args[7]) * 1e6
    window_label <- sprintf("chr%s:%.2f-%.2f Mb", target_chr, start_base/1e6, end_base/1e6)
    out_suffix   <- paste0("_zoom_chr", target_chr)
}

# 3. Filter Region and Execute Keyed Merges
message(paste0("## Filtering genomic profiles across ", window_label))
snp_loc[, chr_clean := gsub("chr", "", as.character(chr), ignore.case = TRUE)]
snp_loc_sub <- snp_loc[chr_clean == target_chr & pos >= start_base & pos <= end_base]

setkey(snp_raw, snpid); setkey(lead_raw, snpid); setkey(snp_loc_sub, snpid)
snp_zoom  <- snp_raw[snp_loc_sub, nomatch = NULL]
lead_zoom <- lead_raw[snp_loc_sub, nomatch = NULL]

if (nrow(snp_zoom) == 0) {
    stop("## ERROR: Zero variants identified inside the target coordinate boundaries.")
}

extreme_hits <- lead_zoom[order(log10p, decreasing = TRUE)][1:min(6, .N)]

# ========================================================================
# PLOT A: DIRECTIONAL COLORING (BETA-BASED)
# ========================================================================
message("## Building Directional Plot (Beta Metrics)...")
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
    scale_x_continuous(expand=c(0.02, 0.02)) + scale_y_continuous(expand=c(0.02, 0.8)) +
    labs(title=paste("Directional Regional Locus Zoom -", basename(out_prefix)),
         subtitle=sprintf("Region: %s | Red (+Beta), Blue (-Beta)", window_label),
         x=paste("Chromosome", target_chr, "Position (Mb)"), y=expression(-log[10](p-value))) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face="bold", size=9)) +
    # FACET STEP: Separates genes into individual clean horizontal rows
    facet_wrap(~geneid, ncol=1, scales="free_y")

if(!is.null(anchor_snp)) {
    p_dir <- p_dir + geom_vline(xintercept=snp_pos/1e6, linetype="dotted", colour="purple", linewidth=0.8)
}

out_file_dir <- paste0(out_prefix, out_suffix, "_directional.png")
png(out_file_dir, width=11, height=6, units="in", res=300)
print(p_dir)
dev.off()

# ========================================================================
# PLOT B: DIAGNOSTIC COLORING BY UNIQUE GENE ID
# ========================================================================
message("## Building Diagnostic Plot (Gene Modulations)...")
snp_bg   <- snp_zoom[fdr_val > fdr_cutoff]
snp_sig  <- snp_zoom[fdr_val <= fdr_cutoff]
lead_sig <- lead_zoom[fdr_val <= fdr_cutoff]

p_gene <- ggplot() +
    geom_point(data=snp_bg, aes(x=pos / 1e6, y=log10p), colour="#999999", size=1.0, alpha=0.3) +
    geom_point(data=snp_sig, aes(x=pos / 1e6, y=log10p, colour=geneid), size=1.4, alpha=0.8) +
    geom_point(data=lead_sig, aes(x=pos / 1e6, y=log10p, colour=geneid), shape=18, size=3.5) +
    geom_hline(yintercept=5, linetype="dashed", colour="grey40", linewidth=0.5) +
    scale_x_continuous(expand=c(0.02, 0.02)) + scale_y_continuous(expand=c(0.02, 0.8)) +
    labs(title=paste("Gene Diagnostic Regional Locus Zoom -", basename(out_prefix)),
         subtitle=sprintf("Region: %s", window_label),
         x=paste("Chromosome", target_chr, "Position (Mb)"), y=expression(-log[10](p-value))) +
    theme_bw() + 
    theme(legend.position="none", # Legend is redundant now since facet headers label the genes
          panel.grid.minor = element_blank(),
          strip.text = element_text(face="bold", size=9)) +
    # FACET STEP: Splits the tracks up cleanly
    facet_wrap(~geneid, ncol=1, scales="free_y")

if(!is.null(anchor_snp)) {
    p_gene p_gene + geom_vline(xintercept=snp_pos/1e6, linetype="dotted", colour="purple", linewidth=0.8)
}

out_file_gene <- paste0(out_prefix, out_suffix, "_by_gene.png")
png(out_file_gene, width=11, height=6.5, units="in", res=300)
print(p_gene)
dev.off()

message(paste("## Track profiles complete!\n## File 1:", out_file_dir, "\n## File 2:", out_file_gene))
