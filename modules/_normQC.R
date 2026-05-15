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
if(!require(Rtsne)) install.packages('Rtsne', repos='http://cran.us.r-project.org')
library(ape)
library(reshape2)
library(Rtsne)

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

# --- REMOVE PROBLEMATIC SUBJECT SAMPLES ---
exclude_samples <- c("CGND-HRA-03028", "CGND-HRA-03029", "CGND-HRA-03030", 
                     "CGND-HRA-03031", "CGND-HRA-03032", "CGND-HRA-03033")
tpm <- tpm[, !colnames(tpm) %in% exclude_samples]

message(paste("# Data dimensions:", nrow(tpm), "genes x", ncol(tpm), "samples"))

# --- Metadata Alignment ---
message("Loading and remapping metadata...")

# Using read.delim with quote="" handles stray quotes in comorbidities or other text fields
covariate <- read.delim(file.path(BASE_DIR, "QTL/covariates.tsv"), 
                        header = TRUE, 
                        sep = "\t", 
                        stringsAsFactors = FALSE, 
                        check.names = FALSE, 
                        quote = "", 
                        fill = TRUE)

# Verification check:
if (nrow(covariate) < 1) {
  stop("Metadata file is empty or could not be parsed correctly.")
}

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

# --- 2. Unrooted Clustering (Fixed & Improved Colors) ---
message("Generating clustering plot...")
sampleDists <- 1 - cor(tpm, method='spearman')
sampleDists[is.na(sampleDists)] <- 1 

hc <- hclust(as.dist(sampleDists), method = "complete")
tree <- as.phylo(hc)

# 1. Extract prefixes (e.g., CGND, SD, GBB)
prefixes <- gsub("^([^-^_]+).*", "\\1", tree$tip.label)
unique_prefixes <- sort(unique(prefixes))

# 2. Use a "Professional" color palette (Darker, non-neon colors)
# We use a color ramp to avoid bright yellows
base_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#000000")
if(length(unique_prefixes) <= length(base_colors)) {
  prefix_colors <- base_colors[1:length(unique_prefixes)]
} else {
  prefix_colors <- colorRampPalette(c("darkblue", "darkred", "darkgreen", "purple4"))(length(unique_prefixes))
}
names(prefix_colors) <- unique_prefixes

# 3. Map prefixes to the tips
tip_colors <- prefix_colors[prefixes]

# 4. Plot
plot(tree, type = "unrooted", cex=.35, lab4ut='axial', 
     tip.color = tip_colors,
     main = paste(TISSUE, "Clustering (Colored by Prefix)"))

# 5. Legend (Fixed: bty='n' removes the box, cex adjusts size)
legend("bottomleft", 
       legend = names(prefix_colors), 
       fill = prefix_colors, 
       cex = 0.7, 
       bty = "n", 
       border = "black")

# --- 3. Sex Gene tSNE Check ---
message("Generating sex-gene tSNE plot...")

SEX_GENE_MAP <- c(
  "ENSG00000229807", # XIST
  "ENSG00000129824", # RPS4Y1
  "ENSG00000280969", # RPS4Y2
  "ENSG00000012817", # KDM5D
  "ENSG00000067048", # DDX3Y
  "ENSG00000114374", # USP9Y
  "ENSG00000183878", # UTY
  "ENSG00000198692"  # EIF1AY
)

clean_rows <- gsub("\\..*", "", rownames(tpm))

sex_idx <- which(clean_rows %in% SEX_GENE_MAP)

if(length(sex_idx) >= 2) {

  sex_tpm <- tpm[sex_idx, , drop=FALSE]

  log_sex_tpm <- log10(sex_tpm + 0.01)

  # transpose so:
  # rows = samples
  # cols = genes
  tsne_input <- t(log_sex_tpm)

  # z-score genes
  tsne_input <- scale(tsne_input)

  # remove zero-variance genes
  keep <- apply(tsne_input, 2, sd, na.rm=TRUE) > 0
  tsne_input <- tsne_input[, keep, drop=FALSE]

  n_samples <- nrow(tsne_input)

  perplexity <- min(30, floor((n_samples - 1) / 3))

  if(perplexity < 2) {
    perplexity <- 2
  }

  set.seed(1234)

  tsne_res <- Rtsne(
    tsne_input,
    dims=2,
    perplexity=perplexity,
    verbose=FALSE,
    max_iter=1000
  )

  plot_cols <- ifelse(sex == "Male", "blue3", "firebrick2")
  plot_cols[is.na(plot_cols)] <- "gray50"

  plot(
    tsne_res$Y[,1],
    tsne_res$Y[,2],
    col=plot_cols,
    pch=19,
    xlab="tSNE 1",
    ylab="tSNE 2",
    main=paste(TISSUE, "Sex Gene tSNE")
  )

  legend(
    "topright",
    legend=c("Male", "Female", "Unknown"),
    col=c("blue3", "firebrick2", "gray50"),
    pch=19,
    bty="n"
  )

  grid(lty="dotted", col="gray80")

  # Optional sample labels for outliers
  male_center <- colMeans(tsne_res$Y[sex == "Male", , drop=FALSE], na.rm=TRUE)
  female_center <- colMeans(tsne_res$Y[sex == "Female", , drop=FALSE], na.rm=TRUE)

  flagged <- c()

  for(i in seq_len(nrow(tsne_res$Y))) {

    if(is.na(sex[i])) {
      next
    }

    pt <- tsne_res$Y[i, ]

    d_male <- sqrt(sum((pt - male_center)^2))
    d_female <- sqrt(sum((pt - female_center)^2))

    if(sex[i] == "Male" && d_female < d_male) {
      flagged <- c(flagged, i)
    }

    if(sex[i] == "Female" && d_male < d_female) {
      flagged <- c(flagged, i)
    }
  }

  if(length(flagged) > 0) {

    text(
      tsne_res$Y[flagged,1],
      tsne_res$Y[flagged,2],
      labels=colnames(tpm)[flagged],
      pos=4,
      cex=0.5
    )

    message(
      paste(
        "Potential sex mismatches:",
        paste(colnames(tpm)[flagged], collapse=", ")
      )
    )
  }

} else {

  message("Warning: Fewer than 2 sex genes found. Skipping tSNE plot.")
}
dev.off()
message(paste("QC Complete. PDF saved to:", outputfile))
