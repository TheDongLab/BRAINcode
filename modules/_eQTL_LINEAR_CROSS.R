###########################################
# Rscript to run interaction eQTL analysis using Matrix eQTL
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Usage: Rscript $PIPELINE_PATH/_eQTL.R snp.txt expression.txt cov.txt output.txt geneloc.txt snploc.txt
# Author: Xianjun Dong (modified by Zachery Wolfe)
# Version: 1.5 (Interaction Filter Edition)
# Date: 7/1/2026
###########################################

.libPaths(c("~/R/libs", .libPaths()))
library('MatrixEQTL')

args<-commandArgs(TRUE)

SNP_file_name=args[1]
expression_file_name=args[2]
covariates_file_name=args[3]
output_file_name=args[4]
gene_location_file_name = args[5]
snp_location_file_name = args[6]

# CRITICAL: Switch model to LINEAR_CROSS for interaction terms
useModel = modelLINEAR_CROSS  

## Settings
pvOutputThreshold = 5e-3;
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 0;   # removed trans eQTL analysis
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

message("## Load gene expression data ...")
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
# CASE-CONTROL MINOR ALLELE CARRIER FILTER (PREVENT PERFECT SEPARATION CRASHES)
# ==============================================================================
# Filter out SNPs that do not have at least 3 minor allele carriers (genotype > 0)
# in BOTH the control and the ALS case group.
# ==============================================================================
message("## Evaluating minor allele carrier distributions across case/control splits...")

full_cov_matrix <- as.matrix(cvrt)
cov_names = rownames(full_cov_matrix)
interaction_idx = which(cov_names == "is_als")

if(length(interaction_idx) == 0) {
    stop("Error: Could not find 'is_als' in the row headers of your covariate file!")
}

is_als_vec <- as.numeric(full_cov_matrix[interaction_idx, ])

# Scan through Slices of the Genotype matrix to flag valid SNPs efficiently
keep_snps_vector <- c()

for(sl in 1:length(snps)) {
    slice_mat <- snps[[sl]]
    
    # Split the matrix into cases and controls up front
    ctrl_mat <- slice_mat[, is_als_vec == 0, drop = FALSE]
    case_mat <- slice_mat[, is_als_vec == 1, drop = FALSE]
    
    # Calculate how many subjects fall into each genotype group (0, 1, 2) per row
    ctrl_0 <- rowSums(abs(ctrl_mat - 0) < 0.1, na.rm = TRUE)
    ctrl_1 <- rowSums(abs(ctrl_mat - 1) < 0.1, na.rm = TRUE)
    ctrl_2 <- rowSums(abs(ctrl_mat - 2) < 0.1, na.rm = TRUE)
    
    case_0 <- rowSums(abs(case_mat - 0) < 0.1, na.rm = TRUE)
    case_1 <- rowSums(abs(case_mat - 1) < 0.1, na.rm = TRUE)
    case_2 <- rowSums(abs(case_mat - 2) < 0.1, na.rm = TRUE)
    
    # Verify that at least 2 distinct genotype buckets have >= 5 individuals
    ctrl_valid_boxes <- (ctrl_0 >= 5) + (ctrl_1 >= 5) + (ctrl_2 >= 5)
    case_valid_boxes <- (case_0 >= 5) + (case_1 >= 5) + (case_2 >= 5)
    
    # Keep only SNPs that meet this baseline layout requirement in both cohorts
    slice_keep <- (ctrl_valid_boxes >= 2) & (case_valid_boxes >= 2)
    keep_snps_vector <- c(keep_snps_vector, slice_keep)
}

# Apply the filtered structural mask directly across Schedulers of SlicedData
snps$RowReorder(which(keep_snps_vector))

message(paste("## Filter Complete: Dropped", sum(!keep_snps_vector), "volatile SNPs."))
message(paste("## Retained", sum(keep_snps_vector), "statistically stable SNPs for interaction engines."))


# ==============================================================================
# INTERACTION SEPARATION STEP
# ==============================================================================
# By default, modelLINEAR_CROSS crosses the SNPs with every row in the 'cvrt' matrix.
# To test ONLY the SNP x is_als interaction, we split them up:
# 1. cvrt: Contains ONLY the main effect variable to cross (is_als)
# 2. cvrt.shared: Contains variables to control for additively (Sex, Age, PCs)
# ==============================================================================

message("## Separating interaction variable (is_als) from main additive covariates...")

# Create a container holding exclusively the interaction row
cvrt_interaction = SlicedData$new()
cvrt_interaction$CreateFromMatrix( full_cov_matrix[interaction_idx, , drop=FALSE] )

# Slice out the rest to act as additive background control variables
additive_indices = setdiff(1:nrow(full_cov_matrix), interaction_idx)
cvrt_additive = SlicedData$new()
cvrt_additive$CreateFromMatrix( full_cov_matrix[additive_indices, , drop=FALSE] )

# Swap standard containers for the engine call
cvrt_to_cross = cvrt_interaction
cvrt_to_adjust = cvrt_additive

# ==============================================================================

message("## Combining background covariates and setting 'is_als' as the last row...")

# 1. Convert SlicedData objects to standard R matrices
mat_adjust <- as.matrix(cvrt_to_adjust)
mat_cross  <- as.matrix(cvrt_to_cross)

# 2. Row-bind them together (background covariates first, interaction term last)
mat_combined <- rbind(mat_adjust, mat_cross)

# 3. Initialize a clean SlicedData object and fill it with the combined matrix
cvrt_combined <- SlicedData$new()
cvrt_combined$CreateFromMatrix(mat_combined)


if(gene_location_file_name != "" && snp_location_file_name != "")
{
    message("## Load gene/SNP location data...")
    snpspos = read.table(snp_location_file_name, header = TRUE, stringsAsFactors = FALSE);
    genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

    message("## Run the cis-/trans-interaction analysis...")
    
    me = Matrix_eQTL_main(
        snps = snps, 
        gene = gene, 
        cvrt = cvrt_combined,          # Combined covariates object
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
    cat('Detected local interaction eQTLs:', '\n');
    show(me$cis$eqtls)
    cat('Detected distant interaction eQTLs:', '\n');
    show(me$trans$eqtls);
} else {
    message("## Run the standard engine interaction analysis...")

    me = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt_combined,         # Combined covariates object
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
    cat('Detected interaction eQTLs:', '\n');
    show(me$all$eqtls);
}

pdf(paste(output_file_name, "pdf", sep="."));
plot(me);
dev.off();
