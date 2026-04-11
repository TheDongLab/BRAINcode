#!/usr/bin/env Rscript
###########################################
# _eQTL_postprocess.R
# Purpose: Post-process raw Matrix eQTL cis output:
#   1. Join SNP chr/pos from snp_location.txt
#   2. Join gene chr/start/end from gene_location.txt
#   3. Flag telomeric SNPs (< 500kb from chromosome start)
#   4. Filter to FDR < 0.05
#   5. Identify lead SNP per gene (lowest p-value)
#   6. STRATEGY: Select top 1000 pairs (500 Extremes + 500 High-Diversity)
###########################################

library(data.table)

args <- commandArgs(TRUE)
cis_file      <- args[1]
snp_loc_file  <- args[2]
gene_loc_file <- args[3]
out_prefix    <- args[4]
fdr_thresh    <- ifelse(is.na(args[5]), 0.05, as.numeric(args[5]))
# We now target 1000 pairs by default for better variety
top_n_total   <- ifelse(is.na(args[6]), 1000, as.integer(args[6]))

TELOMERE_DIST <- 500000   # 500kb

message("## Reading cis-eQTL results...")
eqtl <- fread(cis_file, sep="\t", header=TRUE)
message(sprintf("   %d raw cis-eQTL associations", nrow(eqtl)))

message("## Reading SNP location file...")
snploc <- fread(snp_loc_file, sep="\t", header=TRUE)
setnames(snploc, c("snpid","chr","pos"))
snploc[, pos := as.integer(pos)]

message("## Reading gene location file...")
geneloc <- fread(gene_loc_file, sep="\t", header=TRUE)
setnames(geneloc, c("geneid","gene_chr","gene_start","gene_end"))

message("## Joining SNP coordinates...")
setnames(eqtl, "SNP", "snpid")
eqtl <- merge(eqtl, snploc, by="snpid", all.x=TRUE)

message("## Joining gene coordinates...")
setnames(eqtl, "gene", "geneid")
eqtl <- merge(eqtl, geneloc, by="geneid", all.x=TRUE)

message("## Flagging telomeric SNPs (pos < 500kb from chr start)...")
eqtl[, telomeric_flag := ifelse(!is.na(pos) & pos < TELOMERE_DIST, TRUE, FALSE)]
message(sprintf("   %d associations flagged as telomeric", sum(eqtl$telomeric_flag, na.rm=TRUE)))

message("## Filtering to FDR < ", fdr_thresh, "...")
eqtl_fdr <- eqtl[FDR < fdr_thresh]

message("## Identifying lead SNP per gene...")
eqtl_lead <- eqtl_fdr[order(geneid, `p-value`, pos)]
eqtl_lead <- eqtl_lead[!duplicated(geneid)]

########################################################################
# DUAL-DIRECTION TOP N SELECTION
########################################################################
message("## Selecting top pairs for boxplot (Target: ", top_n_total, ")...")

# Robust SNP matrix lookup
path_parts <- strsplit(out_prefix, "/")[[1]]
tissue_dir_name <- path_parts[length(path_parts) - 4] 
snp_matrix_file <- file.path(dirname(dirname(out_prefix)), paste0("snp_", tissue_dir_name, ".txt"))

if(!file.exists(snp_matrix_file)) {
    parent_dir <- dirname(dirname(out_prefix))
    possible_files <- list.files(parent_dir, pattern = "^snp_.*\\.txt$", full.names = TRUE)
    if(length(possible_files) > 0) snp_matrix_file <- possible_files[1] else stop("SNP matrix not found.")
}

message(paste("## Using SNP matrix for diversity check:", snp_matrix_file))
snp_mat <- fread(snp_matrix_file, header=TRUE)

# CATEGORY 1: The Extremes (Top 500 by P-value regardless of diversity)
# This ensures you get your most significant peaks even if distribution is skewed
eqtl_extremes <- eqtl_fdr[order(`p-value`)][1:min(500, .N), .(geneid, snpid)]

# CATEGORY 2: The High-Diversity Hits (Top 500 with N >= 3 in 2+ groups)
# This captures the 'pretty' plots with lots of samples in each column
candidates <- eqtl_fdr[telomeric_flag == FALSE][order(`p-value`)]
eqtl_diversity <- data.table()

for (i in seq_len(nrow(candidates))) {
    if (nrow(eqtl_diversity) >= 500) break
    
    sid <- candidates$snpid[i]
    geno_vec <- unlist(snp_mat[snpid == sid, -1, with=FALSE])
    counts <- table(factor(geno_vec, levels=0:2))
    
    if (sum(counts >= 3) >= 2) {
        eqtl_diversity <- rbind(eqtl_diversity, candidates[i, .(geneid, snpid)])
    }
}

# Merge Categories
eqtl_top <- unique(rbind(eqtl_extremes, eqtl_diversity))
message(sprintf("   %d total pairs selected (Extremes + High-Diversity)", nrow(eqtl_top)))

########################################################################
# Writing Output Files
########################################################################
message("## Writing output files...")
fwrite(eqtl, file=paste0(out_prefix, ".full_annotated.txt"), sep="\t", quote=FALSE)
fwrite(eqtl_fdr, file=paste0(out_prefix, ".FDR", fdr_thresh, ".txt"), sep="\t", quote=FALSE)
fwrite(eqtl_lead, file=paste0(out_prefix, ".lead_snps.txt"), sep="\t", quote=FALSE)
fwrite(eqtl_top, file=paste0(out_prefix, ".top_for_boxplot.txt"), sep="\t", quote=FALSE, col.names=FALSE)

########################################################################
# Sex-Stratified Meta-Analysis
########################################################################
message("## Checking for meta-analysis opportunity...")
current_sex <- ifelse(grepl("Male", out_prefix), "Male", "Female")
other_sex   <- ifelse(current_sex == "Male", "Female", "Male")
other_sex_file <- gsub(paste0("/", current_sex, "/"), paste0("/", other_sex, "/"), paste0(out_prefix, ".full_annotated.txt"))

if (file.exists(other_sex_file)) {
    message(sprintf("   Found %s results. Generating Meta-Analysis...", other_sex))
    combined_dir <- gsub(paste0("/", current_sex, "/results.*"), "/Combined/results", out_prefix)
    dir.create(combined_dir, showWarnings=FALSE, recursive=TRUE)
    
    eqtl_other <- fread(other_sex_file, sep="\t")
    meta <- merge(eqtl[, .(geneid, snpid, chr, pos, beta_curr = beta, tstat_curr = `t-stat`, pval_curr = `p-value`)],
                  eqtl_other[, .(geneid, snpid, beta_other = beta, tstat_other = `t-stat`, pval_other = `p-value`)],
                  by=c("geneid", "snpid"))
    
    meta[, z_curr  := qnorm(pval_curr / 2, lower.tail=FALSE) * sign(tstat_curr)]
    meta[, z_other := qnorm(pval_other / 2, lower.tail=FALSE) * sign(tstat_other)]
    meta[, z_meta  := (z_curr + z_other) / sqrt(2)]
    meta[, p_meta  := 2 * pnorm(abs(z_meta), lower.tail=FALSE)]
    meta[, FDR_meta := p.adjust(p_meta, method="BH")]
    
    tissue_name <- basename(dirname(dirname(dirname(dirname(out_prefix)))))
    out_meta <- file.path(combined_dir, paste0(tissue_name, "_Combined_Meta_eQTL.txt"))
    fwrite(meta[order(p_meta)], file=out_meta, sep="\t", quote=FALSE)
}
message("## Done.")
