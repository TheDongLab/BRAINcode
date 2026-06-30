#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(ggrepel)

# Hardcoded target window
target_chr <- "17"
start_pos  <- 45465000
end_pos    <- 47055000

message("## Reading data matrix tables...")
snp_raw  <- fread("~/donglab/data/target_ALS/Cerebellum/eQTL/results/Cerebellum_eQTL.cis.txt")
lead_raw <- fread("~/donglab/data/target_ALS/Cerebellum/eQTL/results/Cerebellum_eQTL.lead_snps.txt")
snp_loc  <- fread("~/donglab/data/target_ALS/Cerebellum/eQTL/snp_location.txt", select = c("snpid", "chr", "pos"))

# Standardize headers
setnames(snp_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(lead_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(snp_raw, "gene", "geneid", skip_absent = TRUE)
setnames(lead_raw, "gene", "geneid", skip_absent = TRUE)

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
    labs(title = "Cerebellum regional cis-eQTL Zoom (All Variants)",
         subtitle = sprintf("Region: chr17:45.47-47.05 Mb | Highlighted: %s", top_gene),
         x = "Chromosome 17 Position (Mb)", y = expression(-log[10](p-value))) +
    theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none")

out_file1 <- "~/donglab/data/target_ALS/Cerebellum/eQTL/results/Cerebellum_eQTL.simple_zoom.png"
png(out_file1, width = 10, height = 6, units = "in", res = 300)
print(p1)
dev.off()


# ========================================================================
# PLOT 2: OPTION A (LEAD SNPS ONLY - PRUNED VIEW)
# ========================================================================
message("## Building Plot 2: Option A (Lead variants summary structure)...")

# Group variants by Gene ID and keep only the single strongest sentinel variant per gene
lead_snps_per_gene <- snp_zoom[order(geneid, -log10p), head(.SD, 1), by = geneid]

# Distinct coloring mapping for the top driver vs secondary genes
lead_snps_per_gene[, status := ifelse(geneid == top_gene, "Primary Locus Driver", "Secondary Target Gene")]

p2 <- ggplot(data = lead_snps_per_gene, aes(x = pos / 1e6, y = log10p)) +
    # Reference cutoff line
    geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    # Plot pruned points cleanly categorized by significance tier
    geom_point(aes(colour = status), size = 4.0, alpha = 0.9) +
    scale_colour_manual(values = c("Primary Locus Driver" = "#E41A1C", "Secondary Target Gene" = "#377EB8")) +
    # Repel text label for every distinct gene present in the locus window
    geom_text_repel(aes(label = geneid), size = 3.2, fontface = "bold", 
                    box.padding = 0.4, max.overlaps = 20, cluster_groups = FALSE) +
    scale_x_continuous(limits = c(start_pos/1e6, end_pos/1e6), expand = c(0.02, 0.02)) +
    labs(title = "Cerebellum Lead eQTL Variant Summary (Option A)",
         subtitle = "Pruned view: Showing only the top peak variant per unique gene structure",
         x = "Chromosome 17 Position (Mb)", y = expression(-log[10](p-value))) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), 
          legend.position = "top", 
          legend.title = element_blank())

out_file2 <- "~/donglab/data/target_ALS/Cerebellum/eQTL/results/Cerebellum_eQTL.lead_snps_summary.png"
png(out_file2, width = 10, height = 6, units = "in", res = 300)
print(p2)
dev.off()

message("## Execution completed cleanly!")
message(paste("## Saved Layout 1:", out_file1))
message(paste("## Saved Layout 2:", out_file2))
