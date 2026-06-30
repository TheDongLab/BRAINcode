#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(ggrepel)

# Hardcoded target window
target_chr <- "17"
start_pos  <- 45468000
end_pos    <- 47052000

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

# Find the absolute top signal gene in this window to clear up the color clutter
top_gene <- snp_zoom[order(-log10p)][1, geneid]
message(paste("## Top signal gene identified in this window:", top_gene))

# Split data into: Top Gene vs Background / Other Genes
snp_zoom[, color_cat := ifelse(geneid == top_gene & `p-value` <= 0.05, "Top Gene Locus", "Other/Background")]

# Grab the top lead variants in the region for text labels
extreme_hits <- lead_zoom[order(-log10p)][1:min(3, .N)]

message("## Plotting single-panel region view...")
p <- ggplot() +
    # Background / other genes
    geom_point(data = snp_zoom[color_cat == "Other/Background"], aes(x = pos / 1e6, y = log10p), 
               colour = "#B0B0B0", size = 1.0, alpha = 0.4) +
    # Highlight the primary driver locus safely in one color
    geom_point(data = snp_zoom[color_cat == "Top Gene Locus"], aes(x = pos / 1e6, y = log10p), 
               colour = "#E41A1C", size = 1.4, alpha = 0.8) +
    # Diamonds for lead variants
    geom_point(data = lead_zoom, aes(x = pos / 1e6, y = log10p), 
               shape = 18, size = 3.5, colour = "black") +
    
    geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    
    geom_text_repel(data = extreme_hits, aes(x = pos / 1e6, y = log10p, label = geneid), 
                    size = 3.5, colour = "black", fontface = "bold", box.padding = 0.5) +
    
    scale_x_continuous(limits = c(start_pos/1e6, end_pos/1e6), expand = c(0.01, 0.01)) +
    labs(title = "Cerebellum regional cis-eQTL Zoom",
         subtitle = sprintf("Region: chr17:45.47-47.05 Mb | Highlighted: %s (Top Locus)", top_gene),
         x = "Chromosome 17 Position (Mb)", y = expression(-log[10](p-value))) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), legend.position = "none")

# Save a clean, single 10x6 panel image
out_file <- "~/donglab/data/target_ALS/Cerebellum/eQTL/results/Cerebellum_eQTL.simple_zoom.png"
png(out_file, width = 10, height = 6, units = "in", res = 300)
print(p)
dev.off()
message(paste("## Done! Plot saved to:", out_file))
