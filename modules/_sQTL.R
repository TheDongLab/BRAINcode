###########################################
# Rscript to run sQTL analysis using Matrix eQTL
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Usage: Rscript $PIPELINE_PATH/_sQTL.R snp.txt splicing.txt cov.txt output.txt splicingloc.txt snploc.txt [run_type] [interaction_var]
# Author: Xianjun Dong (modified by Zachery Wolfe)
# Version: 1.2
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

# Read interaction-specific optional arguments if provided
run_type        = if(length(args) >= 7) args[7] else "standard"
interaction_var = if(length(args) >= 8) args[8] else NULL

## Model Configuration
if (run_type == "interaction") {
    message("## Running in INTERACTION mode (modelLINEAR_CROSS)...")
    useModel = modelLINEAR_CROSS
    if (is.null(interaction_var) || interaction_var == "") {
        stop("FATAL ERROR: Interaction mode requested, but no interaction_var parameter was provided.")
    }
} else {
    message("## Running in STANDARD linear mode (modelLINEAR)...")
    useModel = modelLINEAR
}

## Settings
# Only associations significant at this level will be saved
pvOutputThreshold = 5e-3;

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 0;   # removed trans sQTL analysis

# Distance for local intron-SNP pairs
cisDist = 1e6;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

message("## Load genotype data...")

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

message("## Load splicing expression data ...")

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

message("## Load covariates...")

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in slices of 2,000 rows
if(length(covariates_file_name)>0 && file.exists(covariates_file_name)) {
    cvrt$LoadFile(covariates_file_name);
    
    # Process interaction term placement for Matrix eQTL
    if (run_type == "interaction") {
        cov_names = rownames(cvrt)
        if (!interaction_var %in% cov_names) {
            stop(sprintf("FATAL ERROR: Specified interaction variable '%s' not found in covariate file matrix.", interaction_var))
        }
        
        message(sprintf("## Reordering covariate matrix so '%s' is positioned as the last row for interaction testing...", interaction_var))
        
        # Pull data matrix out of SlicedData, reorder rows, and push back
        cvrt_mat = as.matrix(cvrt)
        other_vars = setdiff(rownames(cvrt_mat), interaction_var)
        reordered_mat = cvrt_mat[c(other_vars, interaction_var), , drop = FALSE]
        
        cvrt = SlicedData$new()
        cvrt$CreateFromMatrix(reordered_mat)
    }
} else if (run_type == "interaction") {
    stop("FATAL ERROR: Interaction mode requires a valid covariates file containing the interaction variable.")
}

if(gene_location_file_name != "" && snp_location_file_name!="")
{
    message("## Load splicing/SNP location data...")
    snpspos = read.table(snp_location_file_name, header = TRUE, stringsAsFactors = FALSE);
    genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

    message("## Run the cis-/trans-sQTL analysis...")
    
    me = Matrix_eQTL_main(
    snps = snps, 
    gene = gene, 
    cvrt = cvrt,
    output_file_name     = paste(output_file_name,"trans.txt", sep="."),
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    output_file_name.cis = paste(output_file_name,"cis.txt", sep="."),
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos, 
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
    
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected local sQTLs:', '\n');
    show(me$cis$eqtls)
    cat('Detected distant sQTLs:', '\n');
    show(me$trans$eqtls);
} else {
    message("## Run the sQTL analysis...")

    me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = paste(output_file_name,"txt", sep="."),
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
    
    message("## Getting results ...")

    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected sQTLs:', '\n');
    show(me$all$eqtls);

}

pdf(paste(output_file_name, "pdf", sep="."));
plot(me);
dev.off();
