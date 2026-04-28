###########################################
# Rscript to generate RLE and clustering plots
# Fully Dynamic Path Construction
###########################################

args <- commandArgs(TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript _normQC.R <Tissue_Name>")
}

if(!require(ape)) install.packages('ape', repos='http://cran.us.r-project.org')
if(!require(reshape2)) install.packages('reshape2', repos='http://cran.us.r-project.org')
library(ape)
library(reshape2)

# --- Dynamic Path Construction ---
TISSUE     <- args[1]
BASE_DIR   <- "~/donglab/data/target_ALS"
TISSUE_DIR <- file.path(BASE_DIR, TISSUE, "eQTL")

# Input: joint_TPM_matrix.tsv
TPMfile    <- file.path(TISSUE_DIR, "joint_TPM_matrix.tsv")
# Output: [Tissue]_QC_Plots.pdf
outputfile <- file.path(TISSUE_DIR, paste0(TISSUE, "_QC_Plots.pdf"))

if (!file.exists(TPMfile)) {
  stop(paste("Error: TPM matrix not found at", TPMfile))
}

message(paste("Processing Tissue:", TISSUE))
message(paste("Input File:", TPMfile))
message(paste("Output File:", outputfile))

# Load expression matrix
tpm <- read.table(TPMfile, header=T, check.names=F, row.names=1, sep="\t")
message(paste("# Data dimensions:", nrow(tpm), "genes x", ncol(tpm), "samples"))

# --- Metadata Alignment ---
message("Loading and remapping metadata...")
covariate <- read.table(file.path(BASE_DIR, "QTL/covariates.tsv"), 
                        header=T, sep="\t", stringsAsFactors=F, check.names=F)

TISSUE_REMAP <- list(
  'Motor Cortex Lateral' = 'Motor_Cortex', 'Motor Cortex Medial' = 'Motor_Cortex',
  'Lateral Motor Cortex' = 'Motor_Cortex', 'Medial Motor Cortex' = 'Motor_Cortex',
  'Primary Motor Cortex L' = 'Motor_Cortex', 'Primary Motor Cortex M' = 'Motor_Cortex', 
  'Cortex_Motor_Unspecified' = 'Motor_Cortex', 'Cortex_Motor_BA4' = 'Motor_Cortex', 
  'BA4 Motor Cortex' = 'Motor_Cortex', 'Lateral_motor_cortex' = 'Motor_Cortex',
  'Frontal Cortex' = 'Frontal_Cortex', 'Cerebellum' = 'Cerebellum', 
  'Spinal_Cord_Cervical' = 'Cervical_Spinal_Cord', 'Cervical Spinal Cord' = 'Cervical_Spinal_Cord', 
  'Cervical_spinal_cord' = 'Cervical_Spinal_Cord', 'Spinal_cord_Cervical' = 'Cervical_Spinal_Cord', 
  'Lumbar Spinal Cord' = 'Lumbar_Spinal_Cord', 'Thoracic Spinal Cord' = 'Thoracic_Spinal_Cord', 
  'Spinal_Cord_Lumbosacral' = 'Lumbar_Spinal_Cord', 'Lumbosacral_Spinal_Cord' = 'Lumbar_Spinal_Cord', 
  'Lumbar_spinal_cord' = 'Lumbar_Spinal_Cord'
)

covariate$mapped_tissue <- sapply(covariate$tissue, function(x) {
  if (x %in% names(TISSUE_REMAP)) return(TISSUE_REMAP[[x]])
  return(x)
})

# Filter metadata for current tissue and align with matrix columns
meta_sub <- covariate[covariate$mapped_tissue == TISSUE, ]
meta_idx <- match(colnames(tpm), meta_sub$externalsampleid)

# Extract plot annotations
sex     <- meta_sub$sex[meta_idx]
subject <- meta_sub$externalsubjectid[meta_idx]

pdf(outputfile, width=12, height=10)

# --- 1. RLE Plot ---
message("Generating RLE plot...")
logtpm <- log10(tpm + 1e-4)
rle <- logtpm - apply(logtpm, 1, median)
rle_melt <- melt(as.matrix(rle))
colnames(rle_melt) <- c("Gene", "Sample", "Value")

bymedian <- with(rle_melt, reorder(Sample, Value, IQR))
op <- par(mar=c(8,4,4,2))
boxplot(Value ~ bymedian, data=rle_melt, outline=F, las=2, boxwex=1, 
        col='gray80', cex.axis=0.4, main=paste("RLE:", TISSUE), 
        xlab="", ylab="Relative Log Expression", frame=F)
abline(h=0, col="red", lty=2)

# --- 2. Unrooted Clustering ---
message("Generating clustering plot...")
sampleDists <- 1 - cor(tpm, method='spearman')
hc <- hclust(as.dist(sampleDists), method = "complete")
hc$labels <- colnames(tpm) # Keeping sample IDs for easy outlier identification
tree <- as.phylo(hc)

plot(tree, type = "unrooted", cex=.35, lab4ut='axial', 
     main=paste(TISSUE, "Clustering (Spearman)"))

# --- 3. Gender Marker Check ---
message("Generating gender-match plot...")
clean_rows <- gsub("\\..*", "", rownames(tpm))
idxX <- which(clean_rows == "ENSG00000229807")[1] # XIST
idxY <- which(clean_rows == "ENSG00000129824")[1] # RPS4Y1

if(!is.na(idxX) & !is.na(idxY)) {
  plot_data <- data.frame(
    XIST = as.numeric(logtpm[idxX, ]),
    RPS4Y1 = as.numeric(logtpm[idxY, ]),
    Sex = sex
  )
  plot(plot_data$XIST, plot_data$RPS4Y1, 
       xlab="log10(XIST)", ylab="log10(RPS4Y1)", 
       col=ifelse(plot_data$Sex == "Male", "blue", "red"), 
       pch=19, main=paste(TISSUE, "Gender Check"))
  legend("topright", legend=c("Male", "Female"), col=c("blue", "red"), pch=19)
}

dev.off()
message(paste("QC Complete. PDF saved to:", outputfile))
