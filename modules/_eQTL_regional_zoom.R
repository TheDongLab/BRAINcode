#!/usr/bin/env Rscript
########################################################################
# _eQTL_regional_zoom.R
# Generates high-resolution coordinate-mapped single-chromosome zooms.
# Merges raw p-value metrics with physical locations on the fly.
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

# Set fixed path to reference data
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

# Harmonize column names to prepare for join
setnames(snp_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(lead_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(snp_raw, "gene", "geneid", skip_absent = TRUE)
setnames(lead_raw, "gene", "geneid", skip_absent = TRUE)

# Process metrics
snp_raw[, log10p := -log10(`p-value`)]
snp_raw[, fdr_val := if("FDR" %in% names(snp_raw)) FDR else `p-value`]

lead_raw[, log10p := -log10(`p-value`)]
lead_raw[, fdr_val := if("FDR" %in% names(lead_raw)) FDR else `p-value`]

# 2. Join Coordinates and Filter Region Immediately
message(paste0("## Filtering for Chromosome ", target_chr, "..."))
snp_loc[, chr_clean := gsub("chr", "", as.character(chr), ignore.case = TRUE)]

# Subset map first to speed up the join significantly
snp_loc_sub <- snp_loc[chr_clean == target_chr]

snp_zoom  <- merge(snp_raw, snp_loc_sub, by = "snpid", all.y = FALSE)
lead_zoom <- merge(lead_raw, snp_loc_sub, by = "snpid", all.y = FALSE)

# Filter by window coordinates (converting Mb window inputs to raw bases)
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

# 3. Stratify Layers for Clean ggplot Rendering
snp_bg  <- snp_zoom[fdr_val > fdr_cutoff]
snp_sig <- snp_zoom[fdr_val <= fdr_cutoff]
lead_sig <- lead_zoom[fdr_val <= fdr_cutoff]

# Pick out up to 6 highest peaks to label to prevent messy label crowding
extreme_hits <- lead_zoom[order(log10p, decreasing = TRUE)][1:min(6, .N)]

message("## Building regional visual canvas...")
window_label <- if(!is.null(end_mb)) sprintf("%s:%.2f-%.2f Mb", target_chr, start_mb, end_mb) else sprintf("Chr %s (from %.2f Mb onwards)", target_chr, start_mb)

p3 <- ggplot() +
    # Background variant track
    geom_point(data=snp_bg, aes(x=pos / 1e6, y=log10p), 
               colour="#999999", size=1.0, alpha=0.3) +
    
    # Significant markers highlighted by associated target gene structure
    geom_point(data=snp_sig, aes(x=pos / 1e6, y=log10p, colour=geneid), 
               size=1.4, alpha=0.8) +
    
    # Lead variant diamonds
    geom_point(data=lead_sig, aes(x=pos / 1e6, y=log10p, colour=geneid), 
               shape=18, size=3.5) +
    
    geom_hline(yintercept=-log10(0.05), linetype="dotted", colour="firebrick", linewidth=0.5) +
    geom_hline(yintercept=5, linetype="dashed", colour="grey40", linewidth=0.5) +
    
    geom_text_repel(data=extreme_hits, aes(x=pos / 1e6, y=log10p, label=geneid), 
                    size=3.0, colour="black", fontface="bold", 
                    box.padding = 0.6, max.overlaps = 15) +
    
    scale_x_continuous(expand=c(0.02, 0.02)) +
    scale_y_continuous(expand=c(0.02, 0.8)) +
    
    labs(title=paste("Regional cis-eQTL Target Map -", basename(out_prefix)),
         subtitle=sprintf("Region: %s | Active Loci Groups = %d", window_label, uniqueN(snp_sig$geneid)),
         x=paste("Chromosome", target_chr, "Position (Mb)"), 
         y=expression(-log[10](p-value))) +
    
    theme_bw() + 
    theme(legend.position="bottom", 
          legend.title = element_blank(),
          legend.text = element_text(size = 8.5),
          panel.grid.minor = element_blank())

# Save with localized suffix
out_file <- sprintf("%s.regional_zoom_chr%s.png", out_prefix, target_chr)
png(out_file, width=11, height=6.5, units="in", res=300)
print(p3)
dev.off()

message(paste("## Regional plot successfully exported to:", out_file))
message("## Complete.")
