#!/usr/bin/env Rscript
###########################################
# _eQTL_manhattan.R
# Purpose: Generate two Manhattan plots from annotated cis-eQTL results:
#   1. By SNP position  (standard — where is the variant?)
#   2. By gene position (where is the affected gene?)
#
# Telomeric SNPs are plotted in a distinct colour but not removed.
# Lead SNPs (one per gene) are highlighted with a diamond shape.
#
# Usage:
#   Rscript _eQTL_manhattan.R \
#       annotated_eqtl.txt \
#       lead_snps.txt \
#       output_prefix \
#       [fdr_threshold=0.05]
###########################################

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(scales)
})

args         <- commandArgs(TRUE)
eqtl_file   <- args[1]   # full annotated file from postprocess
lead_file    <- args[2]   # lead SNPs file from postprocess
out_prefix   <- args[3]
fdr_thresh   <- ifelse(is.na(args[4]), 0.05, as.numeric(args[4]))

message("## Reading data...")
eqtl <- fread(eqtl_file, sep="\t", header=TRUE)
lead <- fread(lead_file,  sep="\t", header=TRUE)

# Restrict to autosomes chr1-22 for plotting
eqtl <- eqtl[chr %in% paste0("chr", 1:22)]
lead <- lead[chr %in% paste0("chr", 1:22)]

eqtl[, chr_num := as.integer(sub("chr","", chr))]
lead[, chr_num := as.integer(sub("chr","", chr))]

# FDR significance line
fdr_pval_cutoff <- max(eqtl[FDR < fdr_thresh]$`p-value`, na.rm=TRUE)
sig_line <- -log10(fdr_pval_cutoff)
message(sprintf("   FDR < %.2f corresponds to -log10(p) > %.2f", fdr_thresh, sig_line))

##############################################
# Helper: build cumulative x-axis coordinates
##############################################
build_cumulative <- function(dt, pos_col, chr_col = "chr_num") {
    dt <- copy(dt)
    setorderv(dt, c(chr_col, pos_col))

    # Chromosome max positions for offset calculation
    chr_max <- dt[, .(max_pos = max(get(pos_col), na.rm=TRUE)), by=chr_col]
    setorder(chr_max, chr_num)
    chr_max[, offset := cumsum(shift(max_pos, fill=0))]

    dt <- merge(dt, chr_max[, .(chr_num, offset)], by=chr_col)
    dt[, cum_pos := get(pos_col) + offset]

    # Mid-point of each chromosome for x-axis labels
    axis_df <- dt[, .(centre = mean(range(cum_pos, na.rm=TRUE))), by=chr_col]
    setorder(axis_df, chr_num)

    list(data=dt, axis=axis_df)
}

##############################################
# Shared plot theme
##############################################
manhattan_theme <- function() {
    theme_bw(base_size=13) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "none",
        axis.text.x        = element_text(angle=45, hjust=1, size=9),
        plot.title         = element_text(face="bold", size=13)
    )
}

chr_colours <- c("0"="#555555", "1"="#AAAAAA")   # alternating dark/light

##############################################
# PLOT 1: Manhattan by SNP position
##############################################
message("## Building SNP-position Manhattan plot...")

snp_build  <- build_cumulative(eqtl[!is.na(pos)], "pos")
snp_data   <- snp_build$data
snp_axis   <- snp_build$axis

# Downsample non-significant points for plotting efficiency
set.seed(42)
snp_sig    <- snp_data[`p-value` <= fdr_pval_cutoff]
snp_nonsig <- snp_data[`p-value` >  fdr_pval_cutoff & `p-value` <= 0.001]
snp_nonsig <- snp_nonsig[sample(.N, min(.N, 200000))]
snp_plot   <- rbind(snp_sig, snp_nonsig)
snp_plot[, log10p      := -log10(`p-value`)]
snp_plot[, colour_grp  := as.character(chr_num %% 2)]
snp_plot[, telo_label  := ifelse(telomeric_flag, "telomeric", "normal")]

# Lead SNP positions
lead_snp <- build_cumulative(lead[!is.na(pos)], "pos")$data
lead_snp[, log10p := -log10(`p-value`)]

p1 <- ggplot(snp_plot, aes(x=cum_pos, y=log10p)) +
    geom_point(data=snp_plot[telo_label=="normal"],
               aes(colour=colour_grp), size=0.6, alpha=0.7) +
    geom_point(data=snp_plot[telo_label=="telomeric"],
               colour="#E69F00", size=0.6, alpha=0.9) +
    geom_point(data=lead_snp,
               aes(x=cum_pos, y=log10p),
               shape=18, size=3, colour="red") +
    geom_hline(yintercept=sig_line, linetype="dashed",
               colour="red", linewidth=0.5) +
    scale_colour_manual(values=chr_colours) +
    scale_x_continuous(labels=paste0("chr", snp_axis$chr_num),
                       breaks=snp_axis$centre,
                       expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) +
    labs(title="cis-eQTL Manhattan Plot (by SNP position)",
         x="Chromosome", y=expression(-log[10](p))) +
    annotate("text", x=min(snp_plot$cum_pos), y=sig_line+0.3,
             label=sprintf("FDR < %.2f", fdr_thresh),
             hjust=0, size=3, colour="red") +
    manhattan_theme()

out_snp <- paste0(out_prefix, ".manhattan_by_SNP.pdf")
ggsave(out_snp, p1, width=14, height=5, useDingbats=FALSE)
message(sprintf("   SNP-position Manhattan -> %s", out_snp))

##############################################
# PLOT 2: Manhattan by gene position
##############################################
message("## Building gene-position Manhattan plot...")

gene_data <- eqtl[!is.na(gene_chr) & gene_chr %in% paste0("chr",1:22)]
gene_data[, chr_num   := as.integer(sub("chr","", gene_chr))]
gene_data[, gene_mid  := as.integer((gene_start + gene_end) / 2)]

gene_build <- build_cumulative(gene_data, "gene_mid", "chr_num")
gene_plot  <- gene_build$data
gene_axis  <- gene_build$axis

# Downsample non-significant
gene_sig    <- gene_plot[`p-value` <= fdr_pval_cutoff]
gene_nonsig <- gene_plot[`p-value` >  fdr_pval_cutoff & `p-value` <= 0.001]
gene_nonsig <- gene_nonsig[sample(.N, min(.N, 200000))]
gene_plot2  <- rbind(gene_sig, gene_nonsig)
gene_plot2[, log10p     := -log10(`p-value`)]
gene_plot2[, colour_grp := as.character(chr_num %% 2)]
gene_plot2[, telo_label := ifelse(telomeric_flag, "telomeric", "normal")]

# Lead SNP by gene position
lead_gene <- build_cumulative(
    lead[gene_chr %in% paste0("chr",1:22)][, gene_mid := as.integer((gene_start+gene_end)/2)][, chr_num := as.integer(sub("chr","",gene_chr))],
    "gene_mid", "chr_num")$data
lead_gene[, log10p := -log10(`p-value`)]

p2 <- ggplot(gene_plot2, aes(x=cum_pos, y=log10p)) +
    geom_point(data=gene_plot2[telo_label=="normal"],
               aes(colour=colour_grp), size=0.6, alpha=0.7) +
    geom_point(data=gene_plot2[telo_label=="telomeric"],
               colour="#E69F00", size=0.6, alpha=0.9) +
    geom_point(data=lead_gene,
               aes(x=cum_pos, y=log10p),
               shape=18, size=3, colour="red") +
    geom_hline(yintercept=sig_line, linetype="dashed",
               colour="red", linewidth=0.5) +
    scale_colour_manual(values=chr_colours) +
    scale_x_continuous(labels=paste0("chr", gene_axis$chr_num),
                       breaks=gene_axis$centre,
                       expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01)) +
    labs(title="cis-eQTL Manhattan Plot (by gene position)",
         x="Chromosome", y=expression(-log[10](p))) +
    annotate("text", x=min(gene_plot2$cum_pos), y=sig_line+0.3,
             label=sprintf("FDR < %.2f", fdr_thresh),
             hjust=0, size=3, colour="red") +
    manhattan_theme()

out_gene <- paste0(out_prefix, ".manhattan_by_gene.pdf")
ggsave(out_gene, p2, width=14, height=5, useDingbats=FALSE)
message(sprintf("   Gene-position Manhattan -> %s", out_gene))

message("## Done.")
