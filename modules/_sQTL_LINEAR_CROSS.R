###########################################
# Rscript to run interaction sQTL analysis using Matrix eQTL
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Usage: Rscript $PIPELINE_PATH/_sQTL_LINEAR_CROSS.R snp.txt splicing.txt cov.txt output.txt splicingloc.txt snploc.txt
# Author: Xianjun Dong (modified by Zachery Wolfe)
# Version: 1.7 (Within-Group Invariant PSI & Carrier Filtering)
# Date: 8/4/2026
###########################################

.libPaths(c("~/R/libs", .libPaths()))
library('MatrixEQTL')

args<-commandArgs(TRUE)

SNP_file_name=args[1]
expression_file_name=args[2] # Splicing PSI matrix
covariates_file_name=args[3]
output_file_name=args[4]
gene_location_file_name = args[5] # Splicing location map
snp_location_file_name = args[6]

# CRITICAL: Switch model to LINEAR_CROSS for interaction terms
useModel = modelLINEAR_CROSS  

## Settings
pvOutputThreshold = 5e-3;
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 0;   # removed trans sQTL analysis
cisDist = 1e6;
errorCovariance = numeric();

message("## Load genotype data...")
snps = SlicedData$new();
snps$fileDelimiter = "\t";      
snps$fileOmitCharacters = "NA"; 
snps$fileSkipRows = 1;          
snps$fileSkipColumns = 1;       
snps$fileSliceSize = 2000;      
snps$LoadFile(SNP_file_name);

message("## Load splicing expression data ...")
gene = SlicedData$new();
gene$fileDelimiter = "\t";      
gene$fileOmitCharacters = "NA"; 
gene$fileSkipRows = 1;          
gene$fileSkipColumns = 1;       
gene$fileSliceSize = 2000;      
gene$LoadFile(expression_file_name);

message("## Load covariates...")
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      
cvrt$fileOmitCharacters = "NA"; 
cvrt$fileSkipRows = 1;          
cvrt$fileSkipColumns = 1;       
cvrt$fileSliceSize = 2000;       

if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

# ==============================================================================
# WITHIN-GROUP INVARIANT & BOUNDARY INTRON PRE-FILTERING
# ==============================================================================
message("## Evaluating group-stratified PSI variability across cases and controls...")

full_cov_matrix <- as.matrix(cvrt)
cov_names = rownames(full_cov_matrix)
interaction_idx = which(cov_names == "is_als")

if(length(interaction_idx) == 0) {
    stop("Error: Could not find 'is_als' in the row headers of your covariate file!")
}

is_als_vec <- as.numeric(full_cov_matrix[interaction_idx, ])
ctrl_mask <- (is_als_vec == 0)
case_mask <- (is_als_vec == 1)

keep_introns_vector <- c()

for(sl in 1:length(gene)) {
    slice_mat <- gene[[sl]]

    ctrl_mat <- slice_mat[, ctrl_mask, drop = FALSE]
    case_mat <- slice_mat[, case_mask, drop = FALSE]

    var_ctrl  <- apply(ctrl_mat, 1, var, na.rm = TRUE)
    var_case  <- apply(case_mat, 1, var, na.rm = TRUE)
    mean_ctrl <- apply(ctrl_mat, 1, mean, na.rm = TRUE)
    mean_case <- apply(case_mat, 1, mean, na.rm = TRUE)

    ctrl_non_bnd <- rowSums(ctrl_mat > 0.001 & ctrl_mat < 0.999, na.rm = TRUE) / ncol(ctrl_mat)
    case_non_bnd <- rowSums(case_mat > 0.001 & case_mat < 0.999, na.rm = TRUE) / ncol(case_mat)

    slice_keep <- (!is.na(var_ctrl) & var_ctrl > 1e-4) &
                  (!is.na(var_case) & var_case > 1e-4) &
                  (!is.na(mean_ctrl) & mean_ctrl >= 0.02 & mean_ctrl <= 0.98) &
                  (!is.na(mean_case) & mean_case >= 0.02 & mean_case <= 0.98) &
                  (ctrl_non_bnd >= 0.05) & (case_non_bnd >= 0.05)

    keep_introns_vector <- c(keep_introns_vector, slice_keep)
}

gene$RowReorder(which(keep_introns_vector))

message(paste("## PSI Filter Complete: Dropped", sum(!keep_introns_vector), "group-monomorphic/saturated introns."))
message(paste("## Retained", sum(keep_introns_vector), "variable introns for sQTL interaction modeling."))

# ==============================================================================
# CASE-CONTROL MINOR ALLELE CARRIER FILTER
# ==============================================================================
message("## Evaluating minor allele carrier distributions across case/control splits...")

keep_snps_vector <- c()

for(sl in 1:length(snps)) {
    slice_mat <- snps[[sl]]

    ctrl_mat <- slice_mat[, ctrl_mask, drop = FALSE]
    case_mat <- slice_mat[, case_mask, drop = FALSE]

    ctrl_0 <- rowSums(abs(ctrl_mat - 0) < 0.1, na.rm = TRUE)
    ctrl_1 <- rowSums(abs(ctrl_mat - 1) < 0.1, na.rm = TRUE)
    ctrl_2 <- rowSums(abs(ctrl_mat - 2) < 0.1, na.rm = TRUE)

    case_0 <- rowSums(abs(case_mat - 0) < 0.1, na.rm = TRUE)
    case_1 <- rowSums(abs(case_mat - 1) < 0.1, na.rm = TRUE)
    case_2 <- rowSums(abs(case_mat - 2) < 0.1, na.rm = TRUE)

    ctrl_valid_boxes <- (ctrl_0 >= 5) + (ctrl_1 >= 5) + (ctrl_2 >= 5)
    case_valid_boxes <- (case_0 >= 5) + (case_1 >= 5) + (case_2 >= 5)

    slice_keep <- (ctrl_valid_boxes >= 2) & (case_valid_boxes >= 2)
    keep_snps_vector <- c(keep_snps_vector, slice_keep)
}

snps$RowReorder(which(keep_snps_vector))

message(paste("## SNP Filter Complete: Dropped", sum(!keep_snps_vector), "volatile SNPs."))
message(paste("## Retained", sum(keep_snps_vector), "statistically stable SNPs for interaction engines."))

# ==============================================================================
# INTERACTION SEPARATION STEP
# ==============================================================================
message("## Separating interaction variable (is_als) from main additive covariates...")

cvrt_interaction = SlicedData$new()
cvrt_interaction$CreateFromMatrix(full_cov_matrix[interaction_idx, , drop=FALSE])

additive_indices = setdiff(1:nrow(full_cov_matrix), interaction_idx)
cvrt_additive = SlicedData$new()
cvrt_additive$CreateFromMatrix(full_cov_matrix[additive_indices, , drop=FALSE])

cvrt_to_cross = cvrt_interaction
cvrt_to_adjust = cvrt_additive

message("## Combining background covariates and setting 'is_als' as the last row...")

mat_adjust <- as.matrix(cvrt_to_adjust)
mat_cross <- as.matrix(cvrt_to_cross)
mat_combined <- rbind(mat_adjust, mat_cross)

cvrt_combined <- SlicedData$new()
cvrt_combined$CreateFromMatrix(mat_combined)

# ==============================================================================
# NEAR-PERFECT CIS INTERACTION PRE-FILTER
# ==============================================================================
if(gene_location_file_name != "" && snp_location_file_name != "") {
    message("## Load splicing/SNP location data...")
    snpspos = read.table(snp_location_file_name, header = TRUE, stringsAsFactors = FALSE)
    genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

    message("## Screening cis SNP-intron pairs for near-perfect interaction fits...")

    r2_max <- 1 - 1e-10
    B <- cbind(Intercept = 1, t(mat_combined))
    qr_B <- qr(B)
    Q <- qr.Q(qr_B)[, seq_len(qr_B$rank), drop = FALSE]

    impute_rows <- function(x) {
        missing <- !is.finite(x)
        if(any(missing)) {
            row_means <- rowMeans(x, na.rm = TRUE)
            x[missing] <- row_means[row(x)[missing]]
        }
        x
    }

    bad_snps_vector <- c()

    for(snp_sl in 1:length(snps)) {
        G <- impute_rows(snps[[snp_sl]])
        snp_ids <- rownames(G)
        snp_pos <- snpspos[match(snp_ids, snpspos[, 1]), , drop = FALSE]

        G0 <- G - (G %*% Q) %*% t(Q)
        Z <- sweep(G, 2, is_als_vec, "*")
        Z0 <- Z - (Z %*% Q) %*% t(Q)

        gg <- rowSums(G0^2)
        zg <- rowSums(Z0 * G0)
        z_ss <- rowSums(Z0^2) - zg^2 / gg

        bad_in_slice <- !is.finite(gg) | gg <= 1e-12 |
                        !is.finite(z_ss) | z_ss <= 1e-12

        for(gene_sl in 1:length(gene)) {
            gene_ids <- rownames(gene[[gene_sl]])
            gene_pos <- genepos[match(gene_ids, genepos[, 1]), , drop = FALSE]

            cis_pair <- outer(as.character(snp_pos[, 2]), as.character(gene_pos[, 2]), "==") &
                        outer(as.numeric(snp_pos[, 3]), as.numeric(gene_pos[, 3]) - cisDist, ">=") &
                        outer(as.numeric(snp_pos[, 3]), as.numeric(gene_pos[, 4]) + cisDist, "<=")

            if(!any(cis_pair, na.rm = TRUE)) next

            Y <- impute_rows(gene[[gene_sl]])
            Y0 <- Y - (Y %*% Q) %*% t(Q)

            gy <- G0 %*% t(Y0)
            zy <- Z0 %*% t(Y0)

            numerator <- zy - sweep(gy, 1, zg / gg, "*")
            y_ss <- sweep(gy^2, 1, gg, "/")
            y_ss <- sweep(y_ss, 2, rowSums(Y0^2), function(x, y) y - x)
            denominator <- sweep(y_ss, 1, z_ss, "*")
            partial_r2 <- numerator^2 / denominator

            bad_pair <- cis_pair & (!is.finite(partial_r2) | partial_r2 >= r2_max)
            bad_in_slice <- bad_in_slice | apply(bad_pair, 1, any)
        }

        bad_snps_vector <- c(bad_snps_vector, bad_in_slice)
        message("## Screened SNP slice ", snp_sl, " of ", length(snps))
    }

    snps$RowReorder(which(!bad_snps_vector))

    message(paste("## Interaction R2 Filter Complete: Dropped", sum(bad_snps_vector), "SNPs producing near-perfect cis interaction fits."))
    message(paste("## Retained", sum(!bad_snps_vector), "SNPs for MatrixEQTL."))
}

if(gene_location_file_name != "" && snp_location_file_name != "")
{
    message("## Run the cis-/trans-interaction sQTL analysis...")

    me = Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt_combined,
        output_file_name = paste(output_file_name, "trans.txt", sep="."),
        pvOutputThreshold = pvOutputThreshold_tra,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = paste(output_file_name, "cis.txt", sep="."),
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    );

    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected local interaction sQTLs:', '\n');
    show(me$cis$eqtls)
    cat('Detected distant interaction sQTLs:', '\n');
    show(me$trans$eqtls);
} else {
    message("## Run the standard engine interaction sQTL analysis...")

    me = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt_combined,
        output_file_name = paste(output_file_name, "txt", sep="."),
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    );

    message("## Getting results ...")
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected interaction sQTLs:', '\n');
    show(me$all$eqtls);
}

pdf(paste(output_file_name, "pdf", sep="."));
plot(me);
dev.off();
