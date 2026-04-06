#!/usr/bin/env Rscript
###########################################
# _eQTL_postprocess.R
# Purpose: Post-process raw Matrix eQTL cis output:
#   1. Join SNP chr/pos from snp_location.txt
#   2. Join gene chr/start/end from gene_location.txt
#   3. Flag telomeric SNPs (< 500kb from chromosome start)
#   4. Filter to FDR < 0.05
#   5. Identify lead SNP per gene (lowest p-value)
#   6. Write clean output files for downstream plotting
#
# Usage:
#   Rscript _eQTL_postprocess.R \
#       cis_output.txt \
#       snp_location.txt \
#       gene_location.txt \
#       output_prefix \
#       [fdr_threshold=0.05] \
#       [top_n_boxplot=100]
###########################################

library(data.table)

args <- commandArgs(TRUE)
cis_file      <- args[1]
snp_loc_file  <- args[2]
gene_loc_file <- args[3]
out_prefix    <- args[4]
fdr_thresh    <- ifelse(is.na(args[5]), 0.05, as.numeric(args[5]))
top_n         <- ifelse(is.na(args[6]), 100,  as.integer(args[6]))

TELOMERE_DIST <- 500000   # 500kb

message("## Reading cis-eQTL results...")
eqtl <- fread(cis_file, sep="\t", header=TRUE)
message(sprintf("   %d raw cis-eQTL associations", nrow(eqtl)))

message("## Reading SNP location file...")
snploc <- fread(snp_loc_file, sep="\t", header=TRUE)
# Normalise column names
setnames(snploc, c("snpid","chr","pos"))
snploc[, pos := as.integer(pos)]

message("## Reading gene location file...")
geneloc <- fread(gene_loc_file, sep="\t", header=TRUE)
setnames(geneloc, c("geneid","gene_chr","gene_start","gene_end"))
geneloc[, gene_start := as.integer(gene_start)]
geneloc[, gene_end   := as.integer(gene_end)]

message("## Joining SNP coordinates...")
setnames(eqtl, "SNP", "snpid")
eqtl <- merge(eqtl, snploc, by="snpid", all.x=TRUE)

n_missing_snp <- sum(is.na(eqtl$chr))
if (n_missing_snp > 0)
    message(sprintf("   WARNING: %d SNPs had no location entry", n_missing_snp))

message("## Joining gene coordinates...")
setnames(eqtl, "gene", "geneid")
eqtl <- merge(eqtl, geneloc, by="geneid", all.x=TRUE)

n_missing_gene <- sum(is.na(eqtl$gene_chr))
if (n_missing_gene > 0)
    message(sprintf("   WARNING: %d genes had no location entry", n_missing_gene))

message("## Flagging telomeric SNPs (pos < 500kb from chr start)...")
eqtl[, telomeric_flag := ifelse(!is.na(pos) & pos < TELOMERE_DIST, TRUE, FALSE)]
n_telo <- sum(eqtl$telomeric_flag, na.rm=TRUE)
message(sprintf("   %d associations flagged as telomeric (not removed)", n_telo))

message("## Filtering to FDR < ", fdr_thresh, "...")
eqtl_fdr <- eqtl[FDR < fdr_thresh]
message(sprintf("   %d associations pass FDR < %.2f", nrow(eqtl_fdr), fdr_thresh))

message("## Identifying lead SNP per gene...")
# Lead SNP = lowest p-value per gene; ties broken by position
eqtl_lead <- eqtl_fdr[order(geneid, `p-value`, pos)]
eqtl_lead <- eqtl_lead[!duplicated(geneid)]
message(sprintf("   %d unique genes with a lead SNP", nrow(eqtl_lead)))

message("## Selecting top ", top_n, " gene-SNP pairs for boxplot...")
eqtl_top <- eqtl_fdr[!telomeric_flag == TRUE][order(`p-value`)][seq_len(min(top_n, .N))]
message(sprintf("   %d pairs selected (non-telomeric, lowest p-value)", nrow(eqtl_top)))

message("## Writing output files...")

# Full results with coordinates and flags
out_full <- paste0(out_prefix, ".full_annotated.txt")
fwrite(eqtl, file=out_full, sep="\t", quote=FALSE)
message(sprintf("   Full annotated results -> %s", out_full))

# FDR-filtered results
out_fdr <- paste0(out_prefix, ".FDR", fdr_thresh, ".txt")
fwrite(eqtl_fdr, file=out_fdr, sep="\t", quote=FALSE)
message(sprintf("   FDR-filtered results   -> %s", out_fdr))

# Lead SNP per gene
out_lead <- paste0(out_prefix, ".lead_snps.txt")
fwrite(eqtl_lead, file=out_lead, sep="\t", quote=FALSE)
message(sprintf("   Lead SNPs per gene     -> %s", out_lead))

# Top N for boxplot: just geneid and snpid columns
out_top <- paste0(out_prefix, ".top", top_n, "_for_boxplot.txt")
fwrite(eqtl_top[, .(geneid, snpid)], file=out_top, sep="\t", quote=FALSE, col.names=FALSE)
message(sprintf("   Top %d boxplot pairs   -> %s", top_n, out_top))

message("## Summary:")
message(sprintf("   Raw associations       : %d", nrow(eqtl)))
message(sprintf("   Telomeric (flagged)    : %d", n_telo))
message(sprintf("   FDR < %.2f            : %d", fdr_thresh, nrow(eqtl_fdr)))
message(sprintf("   Unique genes           : %d", nrow(eqtl_lead)))
message(sprintf("   Top %d for boxplot    : %d", top_n, nrow(eqtl_top)))
message("## Done.")
