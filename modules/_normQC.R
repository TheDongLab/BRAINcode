###########################################
# Rscript to generate RLE, clustering, and sex-gene PCA QC plots
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
library(DESeq2)

TISSUE <- args[1]

BASE_DIR <- path.expand(file.path("~/donglab/data", EXPERIMENT))
TISSUE_DIR <- file.path(BASE_DIR, TISSUE, "eQTL")

TPMfile <- file.path(TISSUE_DIR, "joint_TPM_matrix.tsv")
outputfile <- file.path(TISSUE_DIR, paste0(TISSUE, "_QC_Plots.pdf"))
mismatch_file <- file.path(TISSUE_DIR, "potential_sex_mismatches.txt")
genotype_mismatch_file <- file.path(
  TISSUE_DIR,
  "WGS_sex_metadata_mismatches.txt"
)

if(!file.exists(TPMfile)) {
  stop(paste("Error: TPM matrix not found at", TPMfile))
}

message(paste("Experiment:", EXPERIMENT))
message(paste("Processing Tissue:", TISSUE))
message(paste("Input File:", TPMfile))
message(paste("Output File:", outputfile))

tpm <- read.table(
  TPMfile,
  header=TRUE,
  check.names=FALSE,
  row.names=1,
  sep="\t"
)

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

  answerals_metadata_file <- file.path(
    BASE_DIR,
    "aals_dataportal_datatable.csv"
  )

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

  names(answerals_metadata) <- trimws(
    sub("^\ufeff", "", names(answerals_metadata))
  )

  if(nrow(answerals_metadata) < 1 || ncol(answerals_metadata) < 1) {
    stop("AnswerALS metadata file is empty or malformed.")
  }

  id_col <- names(answerals_metadata)[1]
  sex_cols <- which(
    tolower(trimws(names(answerals_metadata))) == "sex"
  )

  if(length(sex_cols) != 1) {
    stop("Expected exactly one Sex column in the AnswerALS metadata.")
  }

  answerals_ids <- trimws(as.character(answerals_metadata[[id_col]]))
  answerals_sex <- tolower(
    trimws(as.character(answerals_metadata[[sex_cols]]))
  )

  answerals_sex <- ifelse(
    answerals_sex == "male", "Male",
    ifelse(answerals_sex == "female", "Female", NA_character_)
  )

  covariate <- data.frame(
    externalsampleid=answerals_ids,
    externalsubjectid=answerals_ids,
    tissue=TISSUE,
    sex=answerals_sex,
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

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(covariate)
)

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

###########################################
# AnswerALS WGS sex-check calls
###########################################

genotype_sex <- setNames(character(0), character(0))

if(EXPERIMENT == "answer_ALS") {

  pfile_prefix <- file.path(
    BASE_DIR,
    "QTL/plink/joint_all_chrs_filtered"
  )

  required_pfile <- paste0(
    pfile_prefix,
    c(".pgen", ".pvar", ".psam")
  )

  missing_pfile <- required_pfile[!file.exists(required_pfile)]
  if(length(missing_pfile) > 0) {
    stop(
      paste(
        "AnswerALS PLINK files are missing:",
        paste(missing_pfile, collapse=", ")
      )
    )
  }

  sexcheck_candidates <- c(
    file.path(
      BASE_DIR,
      "QTL/plink/sexcheck/answerALS_sexcheck.sexcheck"
    ),
    file.path(
      BASE_DIR,
      "QTL/plink/joint_samples_qc.sexcheck"
    ),
    file.path(
      BASE_DIR,
      "QTL/plink/joint_all_chrs_raw.sexcheck"
    )
  )

  sexcheck_candidates <- sexcheck_candidates[
    file.exists(sexcheck_candidates)
  ]

  if(length(sexcheck_candidates) == 0) {
    stop(
      paste(
        "No AnswerALS PLINK sex-check report found. Checked:",
        paste(
          c(
            file.path(
              BASE_DIR,
              "QTL/plink/sexcheck/answerALS_sexcheck.sexcheck"
            ),
            file.path(
              BASE_DIR,
              "QTL/plink/joint_samples_qc.sexcheck"
            ),
            file.path(
              BASE_DIR,
              "QTL/plink/joint_all_chrs_raw.sexcheck"
            )
          ),
          collapse=", "
        )
      )
    )
  }

  sexcheck_file <- sexcheck_candidates[1]
  message(paste("Using WGS sex-check report:", sexcheck_file))

  sexcheck <- read.table(
    sexcheck_file,
    header=TRUE,
    check.names=FALSE,
    comment.char="",
    stringsAsFactors=FALSE
  )

  sexcheck_names <- toupper(sub("^#", "", names(sexcheck)))
  iid_col <- match("IID", sexcheck_names)
  snpsex_col <- match("SNPSEX", sexcheck_names)

  if(is.na(iid_col) || is.na(snpsex_col)) {
    stop(
      "PLINK sex-check report must contain IID and SNPSEX columns."
    )
  }

  snpsex <- trimws(as.character(sexcheck[[snpsex_col]]))
  called_sex <- ifelse(
    snpsex == "1", "Male",
    ifelse(snpsex == "2", "Female", NA_character_)
  )

  genotype_sex <- setNames(
    called_sex,
    trimws(as.character(sexcheck[[iid_col]]))
  )

  message(
    paste(
      "WGS sex-check calls matched to TPM samples:",
      sum(colnames(tpm) %in% names(genotype_sex))
    )
  )
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
  groups <- sub("-.*$", "", tree$tip.label)

  unknown_groups <- setdiff(unique(groups), c("CASE", "CTRL"))
  if(length(unknown_groups) > 0) {
    stop(
      paste(
        "Unexpected AnswerALS sample ID prefixes in clustering plot:",
        paste(unknown_groups, collapse=", ")
      )
    )
  }

  group_colors <- c(CASE="#E41A1C", CTRL="#377EB8")
  tip_colors <- unname(group_colors[groups])

  plot(
    tree,
    type="unrooted",
    cex=.35,
    lab4ut="axial",
    tip.color=tip_colors,
    main=paste(TISSUE, "Clustering (CASE vs CTRL)")
  )

  legend(
    "bottomleft",
    legend=c("CASE", "CTRL"),
    fill=group_colors[c("CASE", "CTRL")],
    cex=0.7,
    bty="n",
    border="black"
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
# 3. Sex-Gene Expression PCA
###########################################

message("Generating sex-gene expression PCA plot...")

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

expression_flagged_samples <- character(0)

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
  expression_flagged_samples <- pca_samples[flagged]

  plot_cols <- ifelse(sex_pca == "Male", "blue3", "firebrick2")
  plot_cols[is.na(plot_cols)] <- "gray50"

  if(EXPERIMENT == "answer_ALS") {
    pca_table <- data.frame(
      Sample=pca_samples,
      MetadataSex=sex_pca,
      WGS_SNPSEX=unname(genotype_sex[pca_samples]),
      PC1=pc1,
      PC2=pc2,
      ExpressionNeighborFlag=ifelse(
        seq_along(pc1) %in% flagged, "YES", "NO"
      )
    )
  } else {
    pca_table <- data.frame(
      Sample=pca_samples,
      Sex=sex_pca,
      PC1=pc1,
      PC2=pc2,
      Flagged=ifelse(seq_along(pc1) %in% flagged, "YES", "NO")
    )
  }

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
    main=paste(TISSUE, "Sex-Gene Expression PCA")
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
      labels=pca_samples[flagged],
      pos=4,
      cex=0.55,
      col=plot_cols[flagged]
    )

    write.table(
      expression_flagged_samples,
      file=mismatch_file,
      quote=FALSE,
      row.names=FALSE,
      col.names=FALSE
    )

    message(
      paste(
        "Potential expression-based sex mismatches:",
        paste(expression_flagged_samples, collapse=", ")
      )
    )

  } else {

    file.create(mismatch_file)
    message("No expression-based sex mismatches detected.")
  }

} else {

  message("Fewer than 2 sex genes found. Skipping expression PCA plot.")
}

###########################################
# AnswerALS: PCA using top 500 HVGs
# Raw counts come from each current donor's normalization.tab.
###########################################

if(EXPERIMENT == "answer_ALS") {
  message("Generating AnswerALS top-500-HVG PCA...")

  processed_dir <- file.path(BASE_DIR, TISSUE, "RNAseq", "Processed")
  count_paths <- file.path(processed_dir, colnames(tpm), "normalization.tab")
  has_counts <- file.exists(count_paths) & file.info(count_paths)$size > 0

  if(any(!has_counts)) {
    message(
      "No normalization.tab for ",
      sum(!has_counts),
      " TPM-matrix samples; they will be omitted from this PCA."
    )
  }

  pca_samples <- colnames(tpm)[has_counts]
  count_paths <- count_paths[has_counts]

  if(length(pca_samples) < 3) {
    stop("Fewer than 3 samples have normalization.tab files for the HVG PCA.")
  }

  count_vectors <- lapply(seq_along(count_paths), function(i) {
    tab <- read.delim(
      count_paths[i],
      header=TRUE,
      sep="\t",
      stringsAsFactors=FALSE,
      check.names=FALSE
    )

    if(!all(c("gene_id", "raw_count") %in% names(tab))) {
      stop(paste("Missing gene_id/raw_count columns in", count_paths[i]))
    }

    gene_id <- sub("\\.[0-9]+$", "", as.character(tab$gene_id))
    raw_count <- suppressWarnings(as.numeric(tab$raw_count))
    keep <- !is.na(gene_id) & gene_id != "" &
      !grepl("^__", gene_id) &
      is.finite(raw_count)

    if(!any(keep)) {
      stop(paste("No usable gene counts in", count_paths[i]))
    }

    summed <- rowsum(
      matrix(raw_count[keep], ncol=1),
      group=gene_id[keep],
      reorder=FALSE
    )

    setNames(as.integer(round(summed[, 1])), rownames(summed))
  })

  all_genes <- sort(unique(unlist(lapply(count_vectors, names))))
  count_matrix <- matrix(
    0L,
    nrow=length(all_genes),
    ncol=length(pca_samples),
    dimnames=list(all_genes, pca_samples)
  )

  for(i in seq_along(count_vectors)) {
    v <- count_vectors[[i]]
    count_matrix[names(v), i] <- v
  }

  # Match the paper's low-expression filter: remove genes averaging <1 raw count/sample.
  keep_genes <- rowMeans(count_matrix) >= 1
  count_matrix <- count_matrix[keep_genes, , drop=FALSE]

  if(nrow(count_matrix) < 500) {
    stop("Fewer than 500 genes remain after the mean raw-count >= 1 filter.")
  }

  col_data <- data.frame(row.names=colnames(count_matrix))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=count_matrix,
    colData=col_data,
    design=~1
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  vst_data <- DESeq2::varianceStabilizingTransformation(dds, blind=TRUE)
  vst_matrix <- SummarizedExperiment::assay(vst_data)

  gene_variance <- apply(vst_matrix, 1, var)
  gene_variance <- gene_variance[is.finite(gene_variance)]
  hvg_genes <- names(sort(gene_variance, decreasing=TRUE))[1:500]

  hvg_pca <- prcomp(
    t(vst_matrix[hvg_genes, , drop=FALSE]),
    center=TRUE,
    scale.=FALSE
  )

  hvg_samples <- rownames(hvg_pca$x)
  hvg_sex <- sex[match(hvg_samples, colnames(tpm))]
  hvg_sex[is.na(hvg_sex)] <- "Unknown"

  point_colors <- ifelse(
    hvg_sex == "Male", "blue3",
    ifelse(hvg_sex == "Female", "firebrick2", "gray50")
  )

  hvg_coords <- data.frame(
    Sample=hvg_samples,
    Sex=hvg_sex,
    PC1=hvg_pca$x[, 1],
    PC2=hvg_pca$x[, 2],
    stringsAsFactors=FALSE
  )

  write.table(
    hvg_coords,
    file=file.path(TISSUE_DIR, paste0(TISSUE, "_top500_HVG_PCA_coordinates.txt")),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

  write.table(
    data.frame(
      gene_id=hvg_genes,
      vst_variance=gene_variance[hvg_genes]
    ),
    file=file.path(TISSUE_DIR, paste0(TISSUE, "_top500_HVGs.txt")),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

  plot(
    hvg_pca$x[, 1],
    hvg_pca$x[, 2],
    col=point_colors,
    pch=19,
    xlab=paste0(
      "PC1 (",
      round(100 * summary(hvg_pca)$importance[2, 1], 1),
      "%)"
    ),
    ylab=paste0(
      "PC2 (",
      round(100 * summary(hvg_pca)$importance[2, 2], 1),
      "%)"
    ),
    main=paste(TISSUE, "Top 500 HVG PCA")
  )

  grid(lty="dotted", col="gray80")

  legend(
    "topright",
    legend=c("Male","Female","Unknown"),
    col=c("blue3","firebrick2","gray50"),
    pch=19,
    bty="n"
  )

  message(
    "Top-500-HVG PCA used ",
    length(hvg_samples),
    " samples and ",
    length(hvg_genes),
    " genes."
  )
}

dev.off()

message(paste("QC Complete. PDF saved to:", outputfile))
message(paste("Expression mismatch file:", mismatch_file))

if(EXPERIMENT == "answer_ALS") {
  message(
    paste(
      "WGS genotype mismatch file:",
      genotype_mismatch_file
    )
  )
}
