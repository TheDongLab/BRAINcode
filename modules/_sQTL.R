###########################################
# Rscript to run sQTL analysis using Matrix eQTL
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Usage: Rscript $PIPELINE_PATH/_sQTL.R snp.txt splicing.txt cov.txt output.prefix splicingloc.txt snploc.txt
# Adapted from eQTL engine v1.2 (Xianjun Dong / Zachery Wolfe)
# Version: 1.0
# Date: 2026-05-20
###########################################

.libPaths(c("~/R/libs", .libPaths()))
library('MatrixEQTL')

args <- commandArgs(TRUE)

SNP_file_name          <- args[1]
expression_file_name   <- args[2] # Splicing PSI matrix
covariates_file_name   <- args[3]
output_file_name       <- args[4] # Output prefix from bash script
gene_location_file_name = args[5] # Splicing location map
snp_location_file_name  = args[6]

useModel = modelLINEAR

## Settings
pvOutputThreshold_cis = 2e-2; # Matches eQTL cutoff exactly
pvOutputThreshold_tra = 0;    # Removed trans sQTL analysis to save memory

# Distance for local intron-SNP pairs (1 Megabase)
cisDist = 1e6;

# Error covariance matrix
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
if(length(covariates_file_name) > 0 && file.exists(covariates_file_name)) {
    cvrt$LoadFile(covariates_file_name);
}

if(gene_location_file_name != "" && snp_location_file_name != "") {
    message("## Load splicing/SNP location data...")
    snpspos = read.table(snp_location_file_name, header = TRUE, stringsAsFactors = FALSE);
    genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

    message("## Run the cis-sQTL analysis...")
    me = Matrix_eQTL_main(
        snps = snps, 
        gene = gene, 
        cvrt = cvrt,
        output_file_name     = paste(output_file_name, "trans.txt", sep="."),
        pvOutputThreshold    = pvOutputThreshold_tra,
        useModel             = useModel, 
        errorCovariance      = errorCovariance, 
        verbose              = TRUE, 
        output_file_name.cis = paste(output_file_name, "cis.txt", sep="."),
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos              = snpspos, 
        genepos              = genepos,
        cisDist              = cisDist,
        pvalue.hist          = FALSE, # Set to FALSE to prevent trans null-pointer crashes
        min.pv.by.genesnp    = FALSE,
        noFDRsaveMemory      = FALSE
    );
    
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected local sQTLs:', '\n');
    show(me$cis$eqtls)
} else {
    stop("FATAL ERROR: Genomic coordinates are mandatory for local cis-sQTL analysis.")
}
