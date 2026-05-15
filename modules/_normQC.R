###########################################
# Rscript to generate RLE, clustering,
# and sex-gene tSNE QC plots
###########################################

args <- commandArgs(TRUE)

if(length(args) < 1) {
  stop("Usage: Rscript _normQC.R <Tissue_Name>")
}

if(!require(ape)) install.packages("ape", repos="http://cran.us.r-project.org")
if(!require(reshape2)) install.packages("reshape2", repos="http://cran.us.r-project.org")
if(!require(Rtsne)) install.packages("Rtsne", repos="http://cran.us.r-project.org")

library(ape)
library(reshape2)
library(Rtsne)

TISSUE <- args[1]

BASE_DIR <- "~/donglab/data/target_ALS"
TISSUE_DIR <- file.path(BASE_DIR, TISSUE, "eQTL")

TPMfile <- file.path(TISSUE_DIR, "joint_TPM_matrix.tsv")
outputfile <- file.path(TISSUE_DIR, paste0(TISSUE, "_QC_Plots.pdf"))
mismatch_file <- file.path(TISSUE_DIR, "potential_sex_mismatches.txt")

if(!file.exists(TPMfile)) {
  stop(paste("Error: TPM matrix not found at", TPMfile))
}

message(paste("Processing Tissue:", TISSUE))
message(paste("Input File:", TPMfile))
message(paste("Output File:", outputfile))

tpm <- read.table(TPMfile, header=TRUE, check.names=FALSE, row.names=1, sep="\t")

exclude_samples <- c(
  "CGND-HRA-03028","CGND-HRA-03029","CGND-HRA-03030",
  "CGND-HRA-03031","CGND-HRA-03032","CGND-HRA-03033"
)

tpm <- tpm[, !colnames(tpm) %in% exclude_samples]

message(paste("# Data dimensions:", nrow(tpm), "genes x", ncol(tpm), "samples"))

message("Loading and remapping metadata...")

covariate <- read.delim(
  file.path(BASE_DIR, "QTL/covariates.tsv"),
  header=TRUE,
  sep="\t",
  stringsAsFactors=FALSE,
  check.names=FALSE,
  quote="",
  fill=TRUE
)

if(nrow(covariate) < 1) {
  stop("Metadata file is empty or malformed.")
}

TISSUE_REMAP <- list(
  "Motor Cortex Lateral"="Motor_Cortex",
  "Motor Cortex Medial"="Motor_Cortex",
  "Lateral Motor Cortex"="Motor_Cortex",
  "Medial Motor Cortex"="Motor_Cortex",
  "Primary Motor Cortex L"="Motor_Cortex",
  "Primary Motor Cortex M"="Motor_Cortex",
  "Cortex_Motor_Unspecified"="Motor_Cortex",
  "Cortex_Motor_BA4"="Motor_Cortex",
  "BA4 Motor Cortex"="Motor_Cortex",
  "Lateral_motor_cortex"="Motor_Cortex",
  "Frontal Cortex"="Frontal_Cortex",
  "Cerebellum"="Cerebellum",
  "Spinal_Cord_Cervical"="Cervical_Spinal_Cord",
  "Cervical Spinal Cord"="Cervical_Spinal_Cord",
  "Cervical_spinal_cord"="Cervical_Spinal_Cord",
  "Spinal_cord_Cervical"="Cervical_Spinal_Cord",
  "Lumbar Spinal Cord"="Lumbar_Spinal_Cord",
  "Thoracic Spinal Cord"="Thoracic_Spinal_Cord",
  "Spinal_Cord_Lumbosacral"="Lumbar_Spinal_Cord",
  "Lumbosacral_Spinal_Cord"="Lumbar_Spinal_Cord",
  "Lumbar_spinal_cord"="Lumbar_Spinal_Cord"
)

covariate$mapped_tissue <- sapply(covariate$tissue, function(x) {
  if(x %in% names(TISSUE_REMAP)) return(TISSUE_REMAP[[x]])
  x
})

meta_sub <- covariate[covariate$mapped_tissue == TISSUE, ]
meta_idx <- match(colnames(tpm), meta_sub$externalsampleid)

sex <- meta_sub$sex[meta_idx]
subject <- meta_sub$externalsubjectid[meta_idx]

pdf(outputfile, width=12, height=10)

###########################################
# 1. RLE Plot
###########################################

message("Generating RLE plot...")

logtpm <- log10(tpm + 1e-4)
rle <- logtpm - apply(logtpm, 1, median)

rle_melt <- melt(as.matrix(rle))
colnames(rle_melt) <- c("Gene","Sample","Value")

bymedian <- with(rle_melt, reorder(Sample, Value, IQR))

par(mar=c(8,4,4,2))

boxplot(
  Value ~ bymedian,
  data=rle_melt,
  outline=FALSE,
  las=2,
  boxwex=1,
  col="gray80",
  cex.axis=0.4,
  main=paste("RLE:", TISSUE),
  xlab="",
  ylab="Relative Log Expression",
  frame=FALSE
)

abline(h=0, col="red", lty=2)

###########################################
# 2. Clustering Plot
###########################################

message("Generating clustering plot...")

sampleDists <- 1 - cor(tpm, method="spearman")
sampleDists[is.na(sampleDists)] <- 1

hc <- hclust(as.dist(sampleDists), method="complete")
tree <- as.phylo(hc)

prefixes <- gsub("^([^-^_]+).*", "\\1", tree$tip.label)
unique_prefixes <- sort(unique(prefixes))

base_colors <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3",
  "#FF7F00","#A65628","#F781BF","#000000"
)

if(length(unique_prefixes) <= length(base_colors)) {
  prefix_colors <- base_colors[1:length(unique_prefixes)]
} else {
  prefix_colors <- colorRampPalette(
    c("darkblue","darkred","darkgreen","purple4")
  )(length(unique_prefixes))
}

names(prefix_colors) <- unique_prefixes
tip_colors <- prefix_colors[prefixes]

plot(
  tree,
  type="unrooted",
  cex=.35,
  lab4ut="axial",
  tip.color=tip_colors,
  main=paste(TISSUE, "Clustering (Colored by Prefix)")
)

legend(
  "bottomleft",
  legend=names(prefix_colors),
  fill=prefix_colors,
  cex=0.7,
  bty="n",
  border="black"
)

###########################################
# 3. Sex Gene tSNE
###########################################

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

  tsne_input <- t(log_sex_tpm)
  tsne_input <- scale(tsne_input)

  keep <- apply(tsne_input, 2, sd, na.rm=TRUE) > 0
  tsne_input <- tsne_input[, keep, drop=FALSE]

  tsne_input <- tsne_input[complete.cases(tsne_input), , drop=FALSE]

  dup_mask <- duplicated(tsne_input)

  if(any(dup_mask)) {
    message(paste("Removing", sum(dup_mask), "duplicate samples before tSNE"))
    tsne_input <- tsne_input[!dup_mask, , drop=FALSE]
  }

  tsne_samples <- rownames(tsne_input)
  sex_tsne <- sex[match(tsne_samples, colnames(tpm))]

  n_samples <- nrow(tsne_input)
  perplexity <- min(30, floor((n_samples - 1) / 3))

  if(perplexity < 2) perplexity <- 2

  set.seed(1234)

  tsne_res <- Rtsne(
    tsne_input,
    dims=2,
    perplexity=perplexity,
    max_iter=2000,
    verbose=FALSE,
    check_duplicates=FALSE
  )

  male_idx <- which(sex_tsne == "Male")
  female_idx <- which(sex_tsne == "Female")

  male_center <- colMeans(tsne_res$Y[male_idx, , drop=FALSE], na.rm=TRUE)
  female_center <- colMeans(tsne_res$Y[female_idx, , drop=FALSE], na.rm=TRUE)

  male_dists <- apply(
    tsne_res$Y[male_idx, , drop=FALSE],
    1,
    function(x) sqrt(sum((x - male_center)^2))
  )

  female_dists <- apply(
    tsne_res$Y[female_idx, , drop=FALSE],
    1,
    function(x) sqrt(sum((x - female_center)^2))
  )

  male_thresh <- median(male_dists, na.rm=TRUE) + 3 * mad(male_dists, na.rm=TRUE)
  female_thresh <- median(female_dists, na.rm=TRUE) + 3 * mad(female_dists, na.rm=TRUE)

  flagged <- c()

  for(i in seq_len(nrow(tsne_res$Y))) {

    if(is.na(sex_tsne[i])) next

    pt <- tsne_res$Y[i, ]

    d_male <- sqrt(sum((pt - male_center)^2))
    d_female <- sqrt(sum((pt - female_center)^2))

    if(
      sex_tsne[i] == "Male" &&
      (
        d_female < d_male ||
        d_male > male_thresh
      )
    ) {
      flagged <- c(flagged, i)
    }

    if(
      sex_tsne[i] == "Female" &&
      (
        d_male < d_female ||
        d_female > female_thresh
      )
    ) {
      flagged <- c(flagged, i)
    }
  }

  flagged <- unique(flagged)

  plot_cols <- ifelse(sex_tsne == "Male", "blue3", "firebrick2")
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

  grid(lty="dotted", col="gray80")

  legend(
    "topright",
    legend=c("Male","Female","Unknown"),
    col=c("blue3","firebrick2","gray50"),
    pch=19,
    bty="n"
  )

  points(
    tsne_res$Y[flagged,1],
    tsne_res$Y[flagged,2],
    pch=1,
    cex=2,
    lwd=2
  )

  if(length(flagged) > 0) {

    flagged_samples <- tsne_samples[flagged]

    text(
      tsne_res$Y[flagged,1],
      tsne_res$Y[flagged,2],
      labels=flagged_samples,
      pos=4,
      cex=0.55
    )

    write.table(
      flagged_samples,
      file=mismatch_file,
      quote=FALSE,
      row.names=FALSE,
      col.names=FALSE
    )

    message(
      paste(
        "Potential sex mismatches:",
        paste(flagged_samples, collapse=", ")
      )
    )

  } else {

    file.create(mismatch_file)
    message("No potential sex mismatches detected.")
  }

} else {

  message("Warning: Fewer than 2 sex genes found. Skipping tSNE plot.")
}

dev.off()

message(paste("QC Complete. PDF saved to:", outputfile))
message(paste("Mismatch file:", mismatch_file))
