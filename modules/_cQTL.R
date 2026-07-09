.libPaths(c("~/R/libs", .libPaths()))
library('MatrixEQTL')

args <- commandArgs(TRUE)

SNP_file_name = args[1]
circ_file_name = args[2]
covariates_file_name = args[3]
output_file_name = args[4]
circ_location_file_name = args[5]
snp_location_file_name = args[6]

useModel = modelLINEAR  # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

## Settings
# Only associations significant at this level will be saved
pvOutputThreshold = 5e-3;

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 0;   # removed trans cQTL analysis

# Distance for local circ-SNP pairs
cisDist = 1e6;

# Error covariance matrix
errorCovariance = numeric();

message("## Load genotype data...")

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

message("## Load circRNA percentage data ...")

circ = SlicedData$new();
circ$fileDelimiter = "\t";      # the TAB character
circ$fileOmitCharacters = "NA"; # denote missing values;
circ$fileSkipRows = 1;          # one row of column labels
circ$fileSkipColumns = 1;       # one column of row labels
circ$fileSliceSize = 2000;      # read file in slices of 2,000 rows
circ$LoadFile(circ_file_name);

message("## Load covariates...")

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;       # read file in slices of 2,000 rows
if(length(covariates_file_name) > 0) {
    cvrt$LoadFile(covariates_file_name);
}

if(circ_location_file_name != "" && snp_location_file_name != "") {
    message("## Load circ/SNP location data...")
    snpspos = read.table(snp_location_file_name, header = TRUE, stringsAsFactors = FALSE);
    circpos = read.table(circ_location_file_name, header = TRUE, stringsAsFactors = FALSE);

    message("## Run the cis-/trans-cQTL analysis...")
    
    me = Matrix_eQTL_main(
        snps = snps, 
        gene = circ, 
        cvrt = cvrt,
        output_file_name     = paste(output_file_name, "trans.txt", sep="."),
        pvOutputThreshold     = pvOutputThreshold_tra,
        useModel = useModel, 
        errorCovariance = errorCovariance, 
        verbose = TRUE, 
        output_file_name.cis = paste(output_file_name, "cis.txt", sep="."),
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos, 
        genepos = circpos,
        cisDist = cisDist,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    );
    
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected local cQTLs:', '\n');
    show(me$cis$eqtls)
    cat('Detected distant cQTLs:', '\n');
    show(me$trans$eqtls);
} else {
    message("## Run the cQTL analysis...")

    me = Matrix_eQTL_engine(
        snps = snps,
        gene = circ,
        cvrt = cvrt,
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
    cat('Detected cQTLs:', '\n');
    show(me$all$eqtls);
}

pdf(paste(output_file_name, "pdf", sep="."));
plot(me);
dev.off();
