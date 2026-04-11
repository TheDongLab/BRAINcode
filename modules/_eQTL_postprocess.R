#!/usr/bin/env Rscript
###########################################
# _eQTL_postprocess.R
###########################################

library(data.table)

args <- commandArgs(TRUE)
cis_file      <- args[1]
snp_loc_file  <- args[2]
gene_loc_file <- args[3]
out_prefix    <- args[4]
fdr_thresh    <- ifelse(is.na(args[5]), 0.05, as.numeric(args[5]))
top_n_total   <- ifelse(is.na(args[6]), 1000, as.integer(args[6]))

TELOMERE_DIST <- 500000 

message("## Reading cis-eQTL results...")
eqtl <- fread(cis_file, sep="\t", header=TRUE)
message("## Reading SNP location file...")
snploc <- fread(snp_loc_file, sep="\t", header=TRUE)
setnames(snploc, c("snpid","chr","pos"))
message("## Reading gene location file...")
geneloc <- fread(gene_loc_file, sep="\t", header=TRUE)
setnames(geneloc, c("geneid","gene_chr","gene_start","gene_end"))

message("## Joining Coordinates...")
setnames(eqtl, "SNP", "snpid")
eqtl <- merge(eqtl, snploc, by="snpid", all.x=TRUE)
setnames(eqtl, "gene", "geneid")
eqtl <- merge(eqtl, geneloc, by="geneid", all.x=TRUE)

eqtl[, telomeric_flag := ifelse(!is.na(pos) & pos < TELOMERE_DIST, TRUE, FALSE)]
eqtl_fdr <- eqtl[FDR < fdr_thresh]

# Lead SNP selection
eqtl_lead <- eqtl_fdr[order(geneid, `p-value`, pos)]
eqtl_lead <- eqtl_lead[!duplicated(geneid)]

########################################################################
# DUAL-DIRECTION TOP N SELECTION
########################################################################
message("## Selecting Top Pairs (Target: ", top_n_total, ")...")

path_parts <- strsplit(out_prefix, "/")[[1]]
tissue_dir_name <- path_parts[length(path_parts) - 4] 
snp_matrix_file <- file.path(dirname(dirname(out_prefix)), paste0("snp_", tissue_dir_name, ".txt"))

if(!file.exists(snp_matrix_file)) {
    parent_dir <- dirname(dirname(out_prefix))
    possible_files <- list.files(parent_dir, pattern = "^snp_.*\\.txt$", full.names = TRUE)
    if(length(possible_files) > 0) snp_matrix_file <- possible_files[1] else stop("SNP matrix not found.")
}

snp_mat <- fread(snp_matrix_file, header=TRUE)

# 1. Extremes (Strongest Signal)
eqtl_extremes <- eqtl_fdr[order(`p-value`)][1:min(500, .N), .(geneid, snpid)]

# 2. Diversity Hits (Ensuring usable boxplots)
candidates <- eqtl_fdr[!(snpid %in% eqtl_extremes$snpid)][telomeric_flag == FALSE][order(`p-value`)]
eqtl_diversity <- data.table()

for (i in seq_len(nrow(candidates))) {
    if (nrow(eqtl_diversity) >= 500) break
    sid <- candidates$snpid[i]
    # Check genotype distribution in SNP matrix
    geno_vec <- unlist(snp_mat[snpid == sid, -1, with=FALSE])
    counts <- table(factor(geno_vec, levels=0:2))
    # Standard: At least 2 groups with at least 3 people
    if (sum(counts >= 3) >= 2) {
        eqtl_diversity <- rbind(eqtl_diversity, candidates[i, .(geneid, snpid)])
    }
}

eqtl_top <- unique(rbind(eqtl_extremes, eqtl_diversity))

fwrite(eqtl, file=paste0(out_prefix, ".full_annotated.txt"), sep="\t", quote=FALSE)
fwrite(eqtl_fdr, file=paste0(out_prefix, ".FDR", fdr_thresh, ".txt"), sep="\t", quote=FALSE)
fwrite(eqtl_lead, file=paste0(out_prefix, ".lead_snps.txt"), sep="\t", quote=FALSE)
fwrite(eqtl_top, file=paste0(out_prefix, ".top_for_boxplot.txt"), sep="\t", quote=FALSE, col.names=FALSE)

# Meta-Analysis Logic
current_sex <- ifelse(grepl("Male", out_prefix), "Male", "Female")
other_sex   <- ifelse(current_sex == "Male", "Female", "Male")
other_sex_file <- gsub(paste0("/", current_sex, "/"), paste0("/", other_sex, "/"), paste0(out_prefix, ".full_annotated.txt"))

if (file.exists(other_sex_file)) {
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
    fwrite(meta[order(p_meta)], file=file.path(combined_dir, paste0(tissue_name, "_Combined_Meta_eQTL.txt")), sep="\t", quote=FALSE)
}
message("## Done.")
