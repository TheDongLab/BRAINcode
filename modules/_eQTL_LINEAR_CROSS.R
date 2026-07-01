###########################################
# Rscript to run interaction eQTL analysis using Matrix eQTL
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Usage: Rscript $PIPELINE_PATH/_eQTL.R snp.txt expression.txt cov.txt output.txt geneloc.txt snploc.txt
# Author: Xianjun Dong (modified by Zachery Wolfe)
# Version: 1.3 (Interaction Edition)
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
# INTERACTION SEPARATION STEP
# ==============================================================================
# By default, modelLINEAR_CROSS crosses the SNPs with every row in the 'cvrt' matrix.
# To test ONLY the SNP x is_als interaction, we split them up:
# 1. cvrt: Contains ONLY the main effect variable to cross (is_als)
# 2. cvrt.shared: Contains variables to control for additively (Sex, Age, PCs)
# ==============================================================================

message("## Separating interaction variable (is_als) from main additive covariates...")

# Find row index of the disease modifier status
cov_names = rownames(cvrt)
interaction_idx = which(cov_names == "is_als")

if(length(interaction_idx) == 0) {
    stop("Error: Could not find 'is_als' in the row headers of your covariate file!")
}

# Create a container holding exclusively the interaction row
cvrt_interaction = cvrt$CreateFromMatrix( cvrt$ToMatrix()[interaction_idx, , drop=FALSE] )

# Slice out the rest to act as additive background control variables
additive_indices = setdiff(1:nrow(cvrt), interaction_idx)
cvrt_additive = cvrt$CreateFromMatrix( cvrt$ToMatrix()[additive_indices, , drop=FALSE] )

# Swap standard containers for the engine call
cvrt_to_cross = cvrt_interaction
cvrt_to_adjust = cvrt_additive

# ==============================================================================

if(gene_location_file_name != "" && snp_location_file_name!="")
{
    message("## Load gene/SNP location data...")
    snpspos = read.table(snp_location_file_name, header = TRUE, stringsAsFactors = FALSE);
    genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

    message("## Run the cis-/trans-interaction analysis...")
    
    me = Matrix_eQTL_main(
        snps = snps, 
        gene = gene, 
        cvrt = cvrt_to_cross,          # Only 'is_als' is crossed with genotype
        cvrt.shared = cvrt_to_adjust,   # Sex, Age, and PCs adjust the model linearly
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
        cvrt = cvrt_to_cross,         # Only 'is_als' is crossed
        cvrt.shared = cvrt_to_adjust,  # Background variations controlled here
        output_file_name = paste(output_file_name,"txt", sep="."),
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
