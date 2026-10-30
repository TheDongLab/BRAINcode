#!/bin/bash
#SBATCH --job-name=sex_gene_heatmap
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/iPSC_Motor_Neuron/eQTL/sex_gene_heatmap_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/iPSC_Motor_Neuron/eQTL/sex_gene_heatmap_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

Rscript - <<'RS'
b <- path.expand("~/donglab/data/answer_ALS")
e <- file.path(b, "iPSC_Motor_Neuron/eQTL")
p <- file.path(b, "iPSC_Motor_Neuron/RNAseq/Processed")
if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("The pheatmap package is required.")
}
g <- c(
    XIST   = "ENSG00000229807",
    JPX    = "ENSG00000225470",
    FTX    = "ENSG00000230590",
    SRY    = "ENSG00000184895",
    RPS4Y1 = "ENSG00000129824",
    RPS4Y2 = "ENSG00000280969",
    DDX3Y  = "ENSG00000067048",
    KDM5D  = "ENSG00000012817",
    UTY    = "ENSG00000183878",
    EIF1AY = "ENSG00000198692",
    USP9Y  = "ENSG00000114374",
    ZFY    = "ENSG00000067646",
    KDM6A  = "ENSG00000147050",
    KDM5C  = "ENSG00000126012",
    DDX3X  = "ENSG00000215301",
    EIF1AX = "ENSG00000173674",
    EIF2S3 = "ENSG00000130741",
    ZFX    = "ENSG00000005889",
    USP9X  = "ENSG00000124486",
    ZRSR2  = "ENSG00000169249"
)
m <- read.csv(
    file.path(b, "aals_dataportal_datatable.csv"),
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM"
)
id_col <- names(m)[1]
sex_col <- which(tolower(trimws(names(m))) == "sex")
if (length(sex_col) != 1) stop("Expected exactly one Sex column.")
sex <- tolower(trimws(as.character(m[[sex_col]])))
sex <- ifelse(sex == "female", "Female",
              ifelse(sex == "male", "Male", "Unknown"))
sex_by_donor <- setNames(sex, trimws(as.character(m[[id_col]])))
mismatch_file <- file.path(e, "potential_sex_mismatches.txt")
if (!file.exists(mismatch_file)) {
    stop(paste("Missing mismatch list:", mismatch_file))
}
mismatches <- unique(trimws(readLines(mismatch_file, warn = FALSE)))
mismatches <- mismatches[nzchar(mismatches)]
if (!length(mismatches)) stop("Mismatch list is empty.")
donor_dirs <- list.dirs(p, recursive = FALSE, full.names = TRUE)
donors <- basename(donor_dirs)
not_found <- setdiff(mismatches, donors)
cat("Loaded", length(mismatches), "unique mismatch IDs.\n")
if (length(not_found)) {
    cat("Mismatch IDs without a donor folder:", paste(not_found, collapse = ", "), "\n")
}
results <- list()
for (d in donor_dirs) {
    file <- file.path(d, "normalization.tab")
    if (!file.exists(file)) next
x <- read.delim(file, check.names = FALSE)
if (!all(c("gene_id", "TPM") %in% names(x))) {
    stop(paste("Missing gene_id or TPM column in", file))
}

gene_ids <- sub("\\.[0-9]+$", "", as.character(x$gene_id))
idx <- match(unname(g), gene_ids)
values <- as.numeric(x$TPM[idx])
names(values) <- names(g)

donor <- basename(d)
donor_sex <- unname(sex_by_donor[donor])
if (is.na(donor_sex) || !length(donor_sex)) donor_sex <- "Unknown"

results[[donor]] <- data.frame(
    donor = donor,
    sex = donor_sex,
    sex_mismatch = ifelse(donor %in% mismatches, "YES", "NO"),
    as.list(values),
    check.names = FALSE
)
}
if (!length(results)) stop("No normalization.tab files were found.")
out <- do.call(rbind, results)
out <- out[order(match(out$sex, c("Female", "Male", "Unknown")), out$donor), ]
if (anyNA(out[, names(g), drop = FALSE])) {
    missing <- names(g)[colSums(is.na(out[, names(g), drop = FALSE])) > 0]
    stop(paste("Missing requested genes in at least one sample:", paste(missing, collapse = ", ")))
}
table_file <- file.path(e, "sex_gene_expression_with_mismatches.tsv")
write.table(out, table_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
heatmap_matrix <- t(as.matrix(out[, names(g), drop = FALSE]))
heatmap_matrix <- log2(heatmap_matrix + 0.01)
heatmap_matrix <- t(scale(t(heatmap_matrix)))
heatmap_matrix[!is.finite(heatmap_matrix)] <- 0
annotation <- data.frame(
    Sex = out$sex,
    SexMismatch = out$sex_mismatch,
    row.names = out$donor
)
female_n <- sum(out$sex == "Female")
known_sex_n <- sum(out$sex %in% c("Female", "Male"))
gaps <- c(female_n, known_sex_n)
gaps <- gaps[gaps > 0 & gaps < nrow(out)]
heatmap_file <- file.path(e, "sex_gene_expression_zscore_heatmap.pdf")
pheatmap::pheatmap(
    heatmap_matrix,
    filename = heatmap_file,
    cluster_cols = FALSE,
    annotation_col = annotation,
    show_colnames = FALSE,
    fontsize_row = 10,
    gaps_col = gaps,
    width = 24,
    height = 8
)
cat(
    "Wrote", nrow(out), "donors; mismatch IDs:", length(mismatches), "\n",
    table_file, "\n",
    heatmap_file, "\n"
)
RS
