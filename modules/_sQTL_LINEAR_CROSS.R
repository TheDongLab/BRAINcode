###########################################
# Rscript to run interaction sQTL analysis using Matrix eQTL
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Usage: Rscript $PIPELINE_PATH/_sQTL_LINEAR_CROSS.R snp.txt splicing.txt cov.txt output.txt splicingloc.txt snploc.txt
# Author: Xianjun Dong (modified by Zachery Wolfe)
# Version: 1.6 (Strategy B - Invariant PSI Filtering)
# Date: 8/3/2026
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

# ==============================================================================
# STRATEGY B: INVARIANT & BOUNDARY INTRON PRE-FILTERING
# ==============================================================================
# Filter out PSI junctions that are invariant or saturated at 0 or 1.
# This eliminates zero-variance residual crashes in modelLINEAR_CROSS.
# ==============================================================================
message("## Evaluating PSI variability and boundary saturation across introns...")

keep_introns_vector <- c()

for(sl in 1:length(gene)) {
    slice_mat <- gene[[sl]]
    
    # Calculate row-wise variance and mean PSI
    row_vars  <- apply(slice_mat, 1, var, na.rm = TRUE)
    row_means <- apply(slice_mat, 1, mean, na.rm = TRUE)
    
    # Calculate proportion of non-boundary samples per intron
    non_boundary_prop <- rowSums(slice_mat > 0.001 & slice_mat < 0.999, na.rm = TRUE) / ncol(slice_mat)
    
    # Keep introns with variance > 1e-3, 0.05 < mean PSI < 0.95, and >= 10% non-boundary samples
    slice_keep <- (!is.na(row_vars) & row_vars > 1e-3) & 
                  (!is.na(row_means) & row_means >= 0.05 & row_means <= 0.95) & 
                  (non_boundary_prop >= 0.10)
    
    keep_introns_vector <- c(keep_introns_vector, slice_keep)
}

# Reorder SlicedData rows to retain only informative introns
gene$RowReorder(which(keep_introns_vector))

message(paste("## PSI Filter Complete: Dropped", sum(!keep_introns_vector), "invariant/saturated introns."))
message(paste("## Retained", sum(keep_introns_vector), "variable introns for sQTL modeling."))

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
# CASE-CONTROL MINOR ALLELE CARRIER FILTER (PREVENT PERFECT SEPARATION CRASHES)
# ==============================================================================
message("## Evaluating minor allele carrier distributions across case/control splits...")

full_cov_matrix <- as.matrix(cvrt)
cov_names = rownames(full_cov_matrix)
interaction_idx = which(cov_names == "is_als")

if(length(interaction_idx) == 0) {
    stop("Error: Could not find 'is_als' in the row headers of your covariate file!")
}

is_als_vec <- as.numeric(full_cov_matrix[interaction_idx, ])

keep_snps_vector <- c()

for(sl in 1:length(snps)) {
    slice_mat <- snps[[sl]]
    
    ctrl_mat <- slice_mat[, is_als_vec == 0, drop = FALSE]
    case_mat <- slice_mat[, is_als_vec == 1, drop = FALSE]
    
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
cvrt_interaction$CreateFromMatrix( full_cov_matrix[interaction_idx, , drop=FALSE] )

additive_indices = setdiff(1:nrow(full_cov_matrix), interaction_idx)
cvrt_additive = SlicedData$new()
cvrt_additive$CreateFromMatrix( full_cov_matrix[additive_indices, , drop=FALSE] )

cvrt_to_cross = cvrt_interaction
cvrt_to_adjust = cvrt_additive

message("## Combining background covariates and setting 'is_als' as the last row...")

mat_adjust <- as.matrix(cvrt_to_adjust)
mat_cross  <- as.matrix(cvrt_to_cross)

mat_combined <- rbind(mat_adjust, mat_cross)

cvrt_combined <- SlicedData$new()
cvrt_combined$CreateFromMatrix(mat_combined)


if(gene_location_file_name != "" && snp_location_file_name != "")
{
    message("## Load splicing/SNP location data...")
    snpspos = read.table(snp_location_file_name, header = TRUE, stringsAsFactors = FALSE);
    genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

    message("## Run the cis-/trans-interaction sQTL analysis...")
    
    me = Matrix_eQTL_main(
        snps = snps, 
        gene = gene, 
        cvrt = cvrt_combined, 
        output_file_name     = paste(output_file_name, "trans.txt", sep="."),
        pvOutputThreshold     = pvOutputThreshold_tra,
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
