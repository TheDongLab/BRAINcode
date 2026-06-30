#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(ggrepel)
library(AnnotationHub)
library(AnnotationDbi)

# ── Parse Command Line Arguments ──────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: No tissue name provided. Usage: Rscript _eQTL_regional_zoom.R <TISSUE>", call. = FALSE)
}
tissue <- args[1]
message(paste("## Running dynamic regional zoom plots for tissue:", tissue))

# Hardcoded target window for Plot 1 & Plot 2
target_chr <- "17"
start_pos  <- 45465000
end_pos    <- 47055000

# Base directory path setup
base_dir <- paste0("~/donglab/data/target_ALS/", tissue, "/eQTL/")
results_dir <- paste0(base_dir, "results/")

message("## Reading data matrix tables...")
snp_raw  <- fread(paste0(results_dir, tissue, "_eQTL.cis.txt"))
lead_raw <- fread(paste0(results_dir, tissue, "_eQTL.lead_snps.txt"))
snp_loc  <- fread(paste0(base_dir, "snp_location.txt"), select = c("snpid", "chr", "pos"))

# Standardize headers
setnames(snp_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(lead_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(snp_raw, "gene", "geneid", skip_absent = TRUE)
setnames(lead_raw, "gene", "geneid", skip_absent = TRUE)
setnames(snp_raw, "source_gene", "geneid", skip_absent = TRUE)
setnames(lead_raw, "source_gene", "geneid", skip_absent = TRUE)

# ── Adaptive Conversion to Common Gene Symbols ────────────────────────
# Check if conversion is actually needed (bypasses DB download if files are already mapped)
needs_conversion <- any(grepl("^ENSG", snp_raw$geneid)) || any(grepl("^ENSG", lead_raw$geneid))

if (needs_conversion) {
  message("## Ensembl IDs detected. Fetching AnnotationHub database for translation...")
  ah <- AnnotationHub()
  human_org <- query(ah, c("OrgDb", "Homo sapiens"))
  orgdb <- ah[[names(human_org)[1]]]
  
  translate_ids <- function(ids) {
    # If the vector is already translated or empty, return as-is
    if (!any(grepl("^ENSG", ids))) return(ids)
    
    # Strip decimals for clean database key-matching
    ids_clean <- gsub("\\..*$", "", as.character(ids))
    
    symbols <- mapIds(orgdb, keys = ids_clean, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
    
    # Fallback to original ID if no mapping symbol exists. NO make.unique() here.
    final_names <- ifelse(is.na(symbols), ids, symbols)
    return(final_names)
  }
  
  message("## Translating gene IDs to common names...")
  if (nrow(snp_raw) > 0) snp_raw[, geneid := translate_ids(geneid)]
  if (nrow(lead_raw) > 0) lead_raw[, geneid := translate_ids(geneid)]
} else {
  message("## Clean common gene names already present in input data matrices. Skipping DB mapping.")
}

# Calculate log p-values
snp_raw[, log10p := -log10(`p-value`)]
lead_raw[, log10p := -log10(`p-value`)]

# Filter location reference map first to save memory overhead
snp_loc[, chr_clean := gsub("chr", "", as.character(chr), ignore.case = TRUE)]
snp_loc_sub <- snp_loc[chr_clean == target_chr & pos >= start_pos & pos <= end_pos]

# Keyed join
setkey(snp_raw, snpid); setkey(lead_raw, snpid); setkey(snp_loc_sub, snpid)
snp_zoom  <- snp_raw[snp_loc_sub, nomatch = NULL]
lead_zoom <- lead_raw[snp_loc_sub, nomatch = NULL]

if (nrow(snp_zoom) == 0) stop("## Zero variants found in this window.")

# Find the absolute top signal gene in this window
top_gene <- snp_zoom[order(-log10p)][1, geneid]
message(paste("## Top signal gene identified in this window:", top_gene))


# ========================================================================
# PLOT 1: THE EXISTING ZOOM PLOT (ALL VARIANTS)
# ========================================================================
message("## Building Plot 1: Full variant zoom...")
snp_zoom[, color_cat := ifelse(geneid == top_gene & `p-value` <= 0.05, "Top Gene Locus", "Other/Background")]
extreme_hits <- lead_zoom[order(-log10p)][1:min(3, .N)]

p1 <- ggplot() +
    geom_point(data = snp_zoom[color_cat == "Other/Background"], aes(x = pos / 1e6, y = log10p), 
               colour = "#B0B0B0", size = 1.0, alpha = 0.4) +
    geom_point(data = snp_zoom[color_cat == "Top Gene Locus"], aes(x = pos / 1e6, y = log10p), 
               colour = "#E41A1C", size = 1.4, alpha = 0.8) +
    geom_point(data = lead_zoom, aes(x = pos / 1e6, y = log10p), 
               shape = 18, size = 3.5, colour = "black") +
    geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text_repel(data = extreme_hits, aes(x = pos / 1e6, y = log10p, label = geneid), 
                    size = 3.5, colour = "black", fontface = "bold", box.padding = 0.5) +
    scale_x_continuous(limits = c(start_pos/1e6, end_pos/1e6), expand = c(0.01, 0.01)) +
    labs(title = sprintf("%s Regional cis-eQTL Zoom (All Variants)", gsub("_", " ", tissue)),
         subtitle = sprintf("Region: chr17:45.47-47.05 Mb | Highlighted: %s", top_gene),
         x = "Chromosome 17 Position (Mb)", y = expression(-log[10](p-value))) +
    theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none")

out_file1 <- paste0(results_dir, tissue, "_eQTL.simple_zoom.png")
png(out_file1, width = 10, height = 6, units = "in", res = 300)
print(p1)
dev.off()


# ========================================================================
# PLOT 2: LEAD SNPS ONLY - PRUNED VIEW
# ========================================================================
message("## Building Plot 2: Option A (Lead variants summary structure)...")

lead_snps_per_gene <- snp_zoom[order(geneid, -log10p), head(.SD, 1), by = geneid]
lead_snps_per_gene[, status := ifelse(geneid == top_gene, "Primary Locus Driver", "Secondary Target Gene")]

p2 <- ggplot(data = lead_snps_per_gene, aes(x = pos / 1e6, y = log10p)) +
    geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_point(aes(colour = status), size = 4.0, alpha = 0.9) +
    scale_colour_manual(values = c("Primary Locus Driver" = "#E41A1C", "Secondary Target Gene" = "#377EB8")) +
    geom_text_repel(aes(label = geneid), size = 3.2, fontface = "bold", 
                    box.padding = 0.4, max.overlaps = 20, cluster_groups = FALSE) +
    scale_x_continuous(limits = c(start_pos/1e6, end_pos/1e6), expand = c(0.02, 0.02)) +
    labs(title = sprintf("%s Lead eQTL Variant Summary", gsub("_", " ", tissue)),
         subtitle = "Pruned view: Showing only the top peak variant per unique gene structure",
         x = "Chromosome 17 Position (Mb)", y = expression(-log[10](p-value))) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), 
          legend.position = "top", 
          legend.title = element_blank())

out_file2 <- paste0(results_dir, tissue, "_eQTL.lead_snps_summary.png")
png(out_file2, width = 10, height = 6, units = "in", res = 300)
print(p2)
dev.off()


# ========================================================================
# PLOT 3: GENE DIAGNOSTIC REGIONAL ZOOM (CENTERED ON rs62056809)
# ========================================================================
message("## Building Plot 3: Gene Diagnostic Zoom centered on rs62056809...")

anchor_snp <- "rs62056809"
anchor_lookup <- snp_loc[snpid == anchor_snp]

if (nrow(anchor_lookup) > 0) {
    anchor_pos <- anchor_lookup[1, pos]
    diag_start <- anchor_pos - 500000
    diag_end   <- anchor_pos + 500000
    
    snp_diag <- snp_raw[snpid %in% snp_loc[chr_clean == target_chr & pos >= diag_start & pos <= diag_end, snpid]]
    setkey(snp_diag, snpid)
    snp_diag <- snp_diag[snp_loc, nomatch = NULL]
    
    sig_genes <- unique(snp_diag[log10p >= 5, geneid])
    snp_diag[, color_group := ifelse(geneid %in% sig_genes & log10p >= 5, geneid, "Background / Non-Sig")]
    
    diag_leads <- lead_raw[snpid %in% snp_diag$snpid & geneid %in% sig_genes]
    diag_leads <- diag_leads[snp_loc, nomatch = NULL, on = "snpid"]
    
    unique_sig_genes <- sort(setdiff(unique(snp_diag$color_group), "Background / Non-Sig"))
    color_palette <- setNames(scales::hue_pal()(length(unique_sig_genes)), unique_sig_genes)
    color_palette["Background / Non-Sig"] <- "#B0B0B0"
    
    p3 <- ggplot() +
        geom_point(data = snp_diag[color_group == "Background / Non-Sig"], aes(x = pos / 1e6, y = log10p),
                   colour = "#B0B0B0", size = 1.0, alpha = 0.3) +
        geom_point(data = snp_diag[color_group != "Background / Non-Sig"], aes(x = pos / 1e6, y = log10p, colour = color_group),
                   size = 1.5, alpha = 0.8) +
        geom_point(data = diag_leads, aes(x = pos / 1e6, y = log10p, colour = geneid),
                   shape = 23, fill = "white", size = 3.0, stroke = 1.5) +
        geom_vline(xintercept = anchor_pos / 1e6, linetype = "dotted", colour = "#984EA3", linewidth = 1.0) +
        geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
        geom_text_repel(data = diag_leads, aes(x = pos / 1e6, y = log10p, label = geneid),
                        size = 3.2, fontface = "bold", box.padding = 0.5, max.overlaps = 15) +
        scale_colour_manual(values = color_palette) +
        scale_x_continuous(limits = c(diag_start / 1e6, diag_end / 1e6), expand = c(0.01, 0.01)) +
        labs(title = sprintf("Gene Diagnostic Regional Zoom - %s eQTL", gsub("_", " ", tissue)),
             subtitle = sprintf("Region: Locus window around %s (chr%s:%.2f-%.2f Mb) | Colored by unique significant target gene structure", 
                                anchor_snp, target_chr, diag_start/1e6, diag_end/1e6),
             x = "Chromosome 17 Position (Mb)", y = expression(-log[10](p-value))) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size = 8)) +
        guides(colour = guide_legend(ncol = 4, byrow = TRUE))
        
    out_file3 <- paste0(results_dir, tissue, "_eQTL.gene_diagnostics.png")
    png(out_file3, width = 11, height = 7, units = "in", res = 300)
    print(p3)
    dev.off()
    message(paste("## Saved Layout 3:", out_file3))
} else {
    message(sprintf("## Warning: Anchor SNP rs62056809 not found in location maps. Skipping Plot 3."))
}

message("## Execution completed cleanly!")
