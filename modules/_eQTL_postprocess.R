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
#   7. Perform Sex Meta-Analysis if both sexes exist
#   8. STRATEGY A: Select top 100 diversity-aware pairs for boxplots
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
# STRATEGY A: Diversity-Aware Top N Selection
########################################################################
message("## Selecting top ", top_n, " gene-SNP pairs for boxplot (N >= 3 check)...")

# Load SNP matrix to check genotype counts
# Path logic: .../eQTL/Male/results/prefix -> .../eQTL/Male/snp_TISSUE.txt
tissue_dir_name <- basename(dirname(dirname(dirname(out_prefix))))
snp_matrix_file <- file.path(dirname(dirname(out_prefix)), paste0("snp_", tissue_dir_name, ".txt"))
snp_mat <- fread(snp_matrix_file, header=TRUE)

# Sort all non-telomeric associations by p-value
candidates <- eqtl_fdr[telomeric_flag == FALSE][order(`p-value`)]
eqtl_top <- data.table()

for (i in seq_len(nrow(candidates))) {
    if (nrow(eqtl_top) >= top_n) break
    
    sid <- candidates$snpid[i]
    # Check diversity: unlist to vector, then table
    geno_vec <- unlist(snp_mat[snpid == sid, -1, with=FALSE])
    counts <- table(factor(geno_vec, levels=0:2))
    
    # Requirement: At least 2 genotype groups have 3 or more samples
    if (sum(counts >= 3) >= 2) {
        eqtl_top <- rbind(eqtl_top, candidates[i, .(geneid, snpid)])
    }
}
message(sprintf("   %d high-diversity pairs selected", nrow(eqtl_top)))

message("## Writing output files...")
out_full <- paste0(out_prefix, ".full_annotated.txt")
fwrite(eqtl, file=out_full, sep="\t", quote=FALSE)

out_fdr <- paste0(out_prefix, ".FDR", fdr_thresh, ".txt")
fwrite(eqtl_fdr, file=out_fdr, sep="\t", quote=FALSE)

out_lead <- paste0(out_prefix, ".lead_snps.txt")
fwrite(eqtl_lead, file=out_lead, sep="\t", quote=FALSE)

out_top <- paste0(out_prefix, ".top", top_n, "_for_boxplot.txt")
fwrite(eqtl_top, file=out_top, sep="\t", quote=FALSE, col.names=FALSE)

########################################################################
# Sex-Stratified Meta-Analysis
########################################################################
message("## Checking for meta-analysis opportunity...")
current_sex <- ifelse(grepl("Male", out_prefix), "Male", "Female")
other_sex   <- ifelse(current_sex == "Male", "Female", "Male")

other_sex_file <- gsub(paste0("/", current_sex, "/"), paste0("/", other_sex, "/"), out_full)
other_sex_file <- gsub(paste0("_", current_sex, "_"), paste0("_", other_sex, "_"), other_sex_file)

if (file.exists(other_sex_file)) {
    message(sprintf("   Found %s results. Generating Combined Meta-Analysis...", other_sex))
    combined_dir <- gsub(paste0("/", current_sex, "/results.*"), "/Combined/results", out_full)
    dir.create(combined_dir, showWarnings=FALSE, recursive=TRUE)
    
    eqtl_other <- fread(other_sex_file, sep="\t")
    meta <- merge(eqtl[, .(geneid, snpid, chr, pos, beta_curr = beta, tstat_curr = `t-stat`, pval_curr = `p-value`)],
                  eqtl_other[, .(geneid, snpid, beta_other = beta, tstat_other = `t-stat`, pval_other = `p-value`)],
                  by=c("geneid", "snpid"))
    
    meta[, z_curr  := qnorm(pval_curr / 2, lower.tail=FALSE) * sign(tstat_curr)]
    meta[, z_other := qnorm(pval_other / 2, lower.tail=FALSE) * sign(tstat_other)]
    meta[, z_meta  := (z_curr + z_other) / sqrt(2)]
    meta[, p_meta  := 2 * pnorm(abs(z_meta), lower.tail=FALSE)]
    meta[, beta_meta := (beta_curr + beta_other) / 2]
    meta[, FDR_meta := p.adjust(p_meta, method="BH")]
    
    tissue_name <- basename(dirname(dirname(dirname(dirname(out_full)))))
    out_meta <- file.path(combined_dir, paste0(tissue_name, "_Combined_Meta_eQTL.txt"))
    fwrite(meta[order(p_meta)], file=out_meta, sep="\t", quote=FALSE)
}
message("## Done.")
