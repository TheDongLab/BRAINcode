###########################################
# Rscript to generate RLE, clustering, and sex-gene tSNE QC plots
###########################################

# Defaults to the legacy TargetALS behavior
# Set EXPERIMENT=answer_ALS in the environment for AnswerALS.
EXPERIMENT <- Sys.getenv("EXPERIMENT", unset="target_ALS")

args <- commandArgs(TRUE)

if(length(args) < 1) {
  stop("Usage: Rscript _normQC.R <Tissue_Name>")
}

if(!EXPERIMENT %in% c("target_ALS", "answer_ALS")) {
  stop("EXPERIMENT must be target_ALS or answer_ALS.")
}

if(!require(ape)) install.packages("ape", repos="http://cran.us.r-project.org")
if(!require(reshape2)) install.packages("reshape2", repos="http://cran.us.r-project.org")
if(!require(Rtsne)) install.packages("Rtsne", repos="http://cran.us.r-project.org")

library(ape)
library(reshape2)
library(Rtsne)

TISSUE <- args[1]

BASE_DIR <- path.expand(file.path("~/donglab/data", EXPERIMENT))
TISSUE_DIR <- file.path(BASE_DIR, TISSUE, "eQTL")

TPMfile <- file.path(TISSUE_DIR, "joint_TPM_matrix.tsv")
outputfile <- file.path(TISSUE_DIR, paste0(TISSUE, "_QC_Plots.pdf"))
mismatch_file <- file.path(TISSUE_DIR, "potential_sex_mismatches.txt")

if(!file.exists(TPMfile)) {
  stop(paste("Error: TPM matrix not found at", TPMfile))
}

message(paste("Experiment:", EXPERIMENT))
message(paste("Processing Tissue:", TISSUE))
message(paste("Input File:", TPMfile))
message(paste("Output File:", outputfile))

tpm <- read.table(TPMfile, header=TRUE, check.names=FALSE, row.names=1, sep="\t")

if(EXPERIMENT == "target_ALS") {
  exclude_samples <- c(
    "CGND-HRA-03028","CGND-HRA-03029","CGND-HRA-03030",
    "CGND-HRA-03031","CGND-HRA-03032","CGND-HRA-03033"
  )

  tpm <- tpm[, !colnames(tpm) %in% exclude_samples, drop=FALSE]
}

message(paste("# Data dimensions:", nrow(tpm), "genes x", ncol(tpm), "samples"))

message("Loading and remapping metadata...")

if(EXPERIMENT == "target_ALS") {

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

} else {

  answerals_metadata_file <- file.path(BASE_DIR, "aals_dataportal_datatable.csv")

  if(!file.exists(answerals_metadata_file)) {
    stop(paste("AnswerALS metadata file not found at", answerals_metadata_file))
  }

  answerals_metadata <- read.csv(
    answerals_metadata_file,
    header=TRUE,
    stringsAsFactors=FALSE,
    check.names=FALSE,
    fileEncoding="UTF-8-BOM"
  )

  names(answerals_metadata) <- trimws(sub("^\ufeff", "", names(answerals_metadata)))

  if(nrow(answerals_metadata) < 1 || ncol(answerals_metadata) < 1) {
    stop("AnswerALS metadata file is empty or malformed.")
  }

  id_col <- names(answerals_metadata)[1]
  sex_cols <- which(tolower(trimws(names(answerals_metadata))) == "sex")
  group_cols <- which(
    gsub("[^[:alnum:]]", "", tolower(names(answerals_metadata))) == "subjectgroup"
  )

  if(length(sex_cols) != 1) {
    stop("Expected exactly one Sex column in the AnswerALS metadata.")
  }

  if(length(group_cols) != 1) {
    stop("Expected exactly one Subject Group column in the AnswerALS metadata.")
  }

  answerals_ids <- trimws(as.character(answerals_metadata[[id_col]]))

  answerals_sex <- tolower(
    trimws(as.character(answerals_metadata[[sex_cols]]))
  )

  answerals_sex <- ifelse(
    answerals_sex == "male", "Male",
    ifelse(answerals_sex == "female", "Female", NA_character_)
  )

  answerals_group <- tolower(
    trimws(as.character(answerals_metadata[[group_cols]]))
  )

  answerals_case_control <- ifelse(
    is.na(answerals_group) | answerals_group == "",
    "Other",
    ifelse(
      answerals_group == "als",
      "Case",
      ifelse(answerals_group == "healthy control", "Control", "Other")
    )
  )

  covariate <- data.frame(
    externalsampleid=answerals_ids,
    externalsubjectid=answerals_ids,
    tissue=TISSUE,
    sex=answerals_sex,
    case_control=answerals_case_control,
    stringsAsFactors=FALSE
  )

  covariate <- covariate[
    !is.na(covariate$externalsampleid) &
      covariate$externalsampleid != "",
    ,
    drop=FALSE
  ]
}

required_metadata_columns <- c(
  "externalsampleid",
  "externalsubjectid",
  "tissue",
  "sex"
)

missing_metadata_columns <- setdiff(required_metadata_columns, names(covariate))
if(length(missing_metadata_columns) > 0) {
  stop(
    paste(
      "Metadata is missing required columns:",
      paste(missing_metadata_columns, collapse=", ")
    )
  )
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

if(TISSUE == "Cortex") {
  meta_sub <- covariate[
    !is.na(covariate$tissue) &
      grepl(
        "(^|[^[:alnum:]])cortex([^[:alnum:]]|$)",
        covariate$tissue,
        ignore.case=TRUE
      ),
  ]
} else if(TISSUE == "Spinal_Cord") {
  meta_sub <- covariate[
    !is.na(covariate$tissue) &
      grepl(
        "(^|[^[:alnum:]])spinal[^[:alnum:]]*cord([^[:alnum:]]|$)",
        covariate$tissue,
        ignore.case=TRUE
      ),
  ]
} else {
  meta_sub <- covariate[covariate$mapped_tissue == TISSUE, ]
}

meta_idx <- match(colnames(tpm), meta_sub$externalsampleid)

unmatched_samples <- colnames(tpm)[is.na(meta_idx)]
if(length(unmatched_samples) > 0) {
  message(
    paste(
      "Warning: TPM samples not matched to metadata:",
      paste(head(unmatched_samples, 20), collapse=", ")
    )
  )
}

sex <- meta_sub$sex[meta_idx]
subject <- meta_sub$externalsubjectid[meta_idx]

if(EXPERIMENT == "answer_ALS") {
  case_control <- meta_sub$case_control[meta_idx]
}

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

if(EXPERIMENT == "answer_ALS") {

  tip_group <- case_control[match(tree$tip.label, colnames(tpm))]
  tip_group[is.na(tip_group)] <- "Other"

  group_colors <- c(
    "Case"="#E41A1C",
    "Control"="#377EB8",
    "Other"="gray60"
  )

  tip_colors <- unname(group_colors[tip_group])

  plot(
    tree,
    type="unrooted",
    cex=.35,
    lab4ut="axial",
    tip.color=tip_colors,
    main=paste(TISSUE, "Clustering (ALS vs Healthy Control)")
  )

  shown_groups <- c("Case", "Control", "Other")
  shown_groups <- shown_groups[shown_groups %in% tip_group]

  legend(
    "bottomleft",
    legend=shown_groups,
    col=group_colors[shown_groups],
    pch=19,
    cex=0.7,
    bty="n"
  )

} else {

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
}

###########################################
# 3. Sex Gene PCA
###########################################

message("Generating sex-gene PCA plot...")

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

pca_table_file <- file.path(
  TISSUE_DIR,
  paste0(TISSUE, "_sex_gene_PCA_coordinates.txt")
)

if(length(sex_idx) >= 2) {

  sex_tpm <- tpm[sex_idx, , drop=FALSE]

  pca_input <- t(log10(sex_tpm + 0.01))
  pca_input <- scale(pca_input)

  keep <- apply(pca_input, 2, sd, na.rm=TRUE) > 0
  pca_input <- pca_input[, keep, drop=FALSE]

  pca_input <- pca_input[complete.cases(pca_input), , drop=FALSE]

  dup_mask <- duplicated(pca_input)

  if(any(dup_mask)) {
    message(paste("Removing", sum(dup_mask), "duplicate samples before PCA"))
    pca_input <- pca_input[!dup_mask, , drop=FALSE]
  }

  pca_samples <- rownames(pca_input)
  sex_pca <- sex[match(pca_samples, colnames(tpm))]

  pca_res <- prcomp(pca_input, center=TRUE, scale.=TRUE)

  pc1 <- pca_res$x[,1]
  pc2 <- pca_res$x[,2]

  coords <- cbind(pc1, pc2)

  flagged <- c()

  for(i in seq_len(nrow(coords))) {

    if(is.na(sex_pca[i])) next

    dists <- sqrt(
      (coords[,1] - coords[i,1])^2 +
      (coords[,2] - coords[i,2])^2
    )

    nn <- order(dists)[2:11]

    nn_sex <- sex_pca[nn]

    male_prop <- mean(nn_sex == "Male", na.rm=TRUE)
    female_prop <- mean(nn_sex == "Female", na.rm=TRUE)

    if(sex_pca[i] == "Female" && male_prop >= 0.8) {
      flagged <- c(flagged, i)
    }

    if(sex_pca[i] == "Male" && female_prop >= 0.8) {
      flagged <- c(flagged, i)
    }
  }

  flagged <- unique(flagged)

  plot_cols <- ifelse(sex_pca == "Male", "blue3", "firebrick2")
  plot_cols[is.na(plot_cols)] <- "gray50"

  pca_table <- data.frame(
    Sample=pca_samples,
    Sex=sex_pca,
    PC1=pc1,
    PC2=pc2,
    Flagged=ifelse(seq_along(pc1) %in% flagged, "YES", "NO")
  )

  write.table(
    pca_table,
    file=pca_table_file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

  plot(
    pc1,
    pc2,
    col=plot_cols,
    pch=19,
    xlab=paste0(
      "PC1 (",
      round(100 * summary(pca_res)$importance[2,1], 1),
      "%)"
    ),
    ylab=paste0(
      "PC2 (",
      round(100 * summary(pca_res)$importance[2,2], 1),
      "%)"
    ),
    main=paste(TISSUE, "Sex Gene PCA")
  )

  grid(lty="dotted", col="gray80")

  legend(
    "topright",
    legend=c("Male","Female","Unknown"),
    col=c("blue3","firebrick2","gray50"),
    pch=19,
    bty="n"
  )

  if(length(flagged) > 0) {

    flagged_samples <- pca_samples[flagged]

    points(
      pc1[flagged],
      pc2[flagged],
      pch=1,
      cex=2,
      lwd=2
    )

    text(
      pc1[flagged],
      pc2[flagged],
      labels=flagged_samples,
      pos=4,
      cex=0.55,
      col=plot_cols[flagged]
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

  message("Warning: Fewer than 2 sex genes found. Skipping PCA plot.")
}

dev.off()

message(paste("QC Complete. PDF saved to:", outputfile))
message(paste("Mismatch file:", mismatch_file))
