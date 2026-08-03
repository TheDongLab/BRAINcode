###########################################
# Rscript to run interaction sQTL analysis using Matrix eQTL
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Usage: Rscript $PIPELINE_PATH/_sQTL_LINEAR_CROSS.R snp.txt splicing.txt cov.txt output.txt splicingloc.txt snploc.txt
# Author: Xianjun Dong (modified by Zachery Wolfe)
# Version: 1.5 (Interaction Filter Edition - Splicing)
# Date: 7/29/2026
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
# CASE-CONTROL INTERACTION STABILITY FILTER (PREVENT CEILING / SEPARATION ARTIFACTS)
# ==============================================================================
# Verifies that:
#   1. Both Controls and Cases have >= 5 individuals in at least 2 genotype buckets.
#   2. Both Controls and Cases have non-zero phenotype variance (var > 1e-4) to prevent flatline baseline slices (e.g. all PSI = 0) from driving SE -> 0.
# ==============================================================================
message("## Evaluating genotype distributions and per-group phenotype variance across case/control splits...")

full_cov_matrix <- as.matrix(cvrt)
cov_names <- rownames(full_cov_matrix)
interaction_idx <- which(cov_names == "is_als")

if (length(interaction_idx) == 0) {
    stop("Error: Could not find 'is_als' in the row headers of your covariate file!")
}

is_als_vec <- as.numeric(full_cov_matrix[interaction_idx, ])
ctrl_idx <- which(is_als_vec == 0)
case_idx <- which(is_als_vec == 1)

# ------------------------------------------------------------------------------
# STEP 1: Aggressive Phenotype Filtering
# ------------------------------------------------------------------------------
message("## Filtering splicing phenotype matrix for per-group variance & continuous spread...")

keep_gene_vector <- c()

for (gl in 1:length(gene)) {
    gene_mat <- gene[[gl]]
    
    ctrl_psi <- gene_mat[, ctrl_idx, drop = FALSE]
    case_psi <- gene_mat[, case_idx, drop = FALSE]
    
    # Reject introns where > 80% of samples in EITHER cohort are locked at 0 or 1
    ctrl_is_saturated <- apply(ctrl_psi, 1, function(x) mean(x <= 0.01 | x >= 0.99, na.rm=TRUE) > 0.80)
    case_is_saturated <- apply(case_psi, 1, function(x) mean(x <= 0.01 | x >= 0.99, na.rm=TRUE) > 0.80)

    # Require non-negligible standard deviation (SD > 0.02) in both groups
    ctrl_sd <- apply(ctrl_psi, 1, sd, na.rm = TRUE) > 0.02
    case_sd <- apply(case_psi, 1, sd, na.rm = TRUE) > 0.02

    slice_gene_keep <- !ctrl_is_saturated & !case_is_saturated & ctrl_sd & case_sd
    
    keep_gene_vector <- c(keep_gene_vector, slice_gene_keep)
}

gene$RowReorder(which(keep_gene_vector))

message(paste("## Phenotype Filter Complete: Dropped", sum(!keep_gene_vector), "boundary-saturated/invariant introns."))
message(paste("## Retained", sum(keep_gene_vector), "introns with robust within-group PSI distributions."))

# ------------------------------------------------------------------------------
# STEP 2: Filter SNP Matrix for Per-Group Genotype Spread
# ------------------------------------------------------------------------------
keep_snps_vector <- c()

for (sl in 1:length(snps)) {
    slice_mat <- snps[[sl]]
    
    # Split the genotype matrix into cases and controls up front
    ctrl_mat <- slice_mat[, ctrl_idx, drop = FALSE]
    case_mat <- slice_mat[, case_idx, drop = FALSE]
    
    # Count subjects falling into each genotype bucket (0, 1, 2)
    ctrl_0 <- rowSums(abs(ctrl_mat - 0) < 0.1, na.rm = TRUE)
    ctrl_1 <- rowSums(abs(ctrl_mat - 1) < 0.1, na.rm = TRUE)
    ctrl_2 <- rowSums(abs(ctrl_mat - 2) < 0.1, na.rm = TRUE)
    
    case_0 <- rowSums(abs(case_mat - 0) < 0.1, na.rm = TRUE)
    case_1 <- rowSums(abs(case_mat - 1) < 0.1, na.rm = TRUE)
    case_2 <- rowSums(abs(case_mat - 2) < 0.1, na.rm = TRUE)
    
    # Require at least 2 distinct genotype buckets with >= 5 subjects in both groups
    ctrl_valid_boxes <- (ctrl_0 >= 5) + (ctrl_1 >= 5) + (ctrl_2 >= 5)
    case_valid_boxes <- (case_0 >= 5) + (case_1 >= 5) + (case_2 >= 5)
    
    # Extra safety: Ensure non-zero genotype variance per subgroup
    ctrl_snp_var <- apply(ctrl_mat, 1, var, na.rm = TRUE)
    case_snp_var <- apply(case_mat, 1, var, na.rm = TRUE)
    
    slice_snp_keep <- (ctrl_valid_boxes >= 2) & (case_valid_boxes >= 2) & 
                       (!is.na(ctrl_snp_var) & ctrl_snp_var > 1e-4) & 
                       (!is.na(case_snp_var) & case_snp_var > 1e-4)
                       
    keep_snps_vector <- c(keep_snps_vector, slice_snp_keep)
}

snps$RowReorder(which(keep_snps_vector))

message(paste("## Genotype Filter Complete: Dropped", sum(!keep_snps_vector), "volatile/unbalanced SNPs."))
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
    message("## Load splicing/SNP location data...")
    snpspos = read.table(snp_location_file_name, header = TRUE, stringsAsFactors = FALSE);
    genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

    message("## Run the cis-/trans-interaction sQTL analysis...")
    
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
    cat('Detected local interaction sQTLs:', '\n');
    show(me$cis$eqtls)
    cat('Detected distant interaction sQTLs:', '\n');
    show(me$trans$eqtls);
} else {
    message("## Run the standard engine interaction sQTL analysis...")

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
    cat('Detected interaction sQTLs:', '\n');
    show(me$all$eqtls);
}

pdf(paste(output_file_name, "pdf", sep="."));
plot(me);
dev.off();
