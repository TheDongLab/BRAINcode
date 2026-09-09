#!/usr/bin/env Rscript

# Static genotype-stratified LeafViz-style plot for sQTL events.
# No installed leafcutter R package is required.

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  out <- list()
  i <- 1
  while (i <= length(x)) {
    key <- x[i]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(x)) stop("Missing value for ", key)
    out[[sub("^--", "", key)]] <- x[i + 1]
    i <- i + 2
  }
  out
}

a <- parse_args(args)

required <- c(
  "tissue", "tissue-regex", "anchor", "gene",
  "variant-chr", "variant-pos", "variant-id", "variant-ref", "variant-alt",
  "genotypes", "metadata", "data-root", "gtf", "outdir", "prefix"
)

missing_args <- required[!required %in% names(a)]
if (length(missing_args)) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("R package 'ggplot2' is required.")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("R package 'gridExtra' is required.")
  }
})

library(ggplot2)
library(gridExtra)
library(grid)

dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)

cat("Loading metadata...\n")
meta <- read.csv(a$metadata, stringsAsFactors = FALSE, check.names = FALSE)

need_meta <- c("externalsampleid", "externalsubjectid", "tissue")
if (!all(need_meta %in% colnames(meta))) {
  stop(
    "Metadata must contain columns: ",
    paste(need_meta, collapse = ", ")
  )
}

meta <- meta[
  !is.na(meta$tissue) &
    grepl(a[["tissue-regex"]], meta$tissue, ignore.case = TRUE, perl = TRUE),
  ,
  drop = FALSE
]

if (!nrow(meta)) {
  stop("No metadata rows matched tissue pattern for ", a$tissue)
}

# File names/directories use underscores rather than dashes.
meta$sample_dir <- gsub("-", "_", meta$externalsampleid, fixed = TRUE)

cat("Loading VCF genotypes...\n")
gt <- read.delim(a$genotypes, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("externalsubjectid", "GT") %in% colnames(gt))) {
  stop("Genotype table must contain externalsubjectid and GT.")
}

normalize_gt <- function(x) {
  x <- gsub("\\|", "/", x)
  out <- rep(NA_character_, length(x))
  out[x == "0/0"] <- "Ref/Ref"
  out[x %in% c("0/1", "1/0")] <- "Het"
  out[x == "1/1"] <- "Hom Alt"
  out
}

gt$group <- normalize_gt(gt$GT)

meta <- merge(
  meta,
  gt[, c("externalsubjectid", "GT", "group")],
  by = "externalsubjectid",
  all.x = TRUE
)

meta <- meta[!is.na(meta$group), , drop = FALSE]

if (!nrow(meta)) {
  stop("No tissue-matched RNA-seq samples had usable 0/0, 0/1, or 1/1 genotypes.")
}

group_levels <- c("Ref/Ref", "Het", "Hom Alt")
meta$group <- factor(meta$group, levels = group_levels)

# Build an index of all per-sample PSI files across tissue folders.
cat("Indexing LeafCutter PSI files...\n")
psi_files <- Sys.glob(
  file.path(
    a[["data-root"]],
    "*",
    "RNAseq",
    "Processed",
    "*",
    "leafcutter",
    "psi",
    "*.leafcutter.PSI.tsv"
  )
)

if (!length(psi_files)) {
  stop("No LeafCutter PSI files found beneath ", a[["data-root"]])
}

# Prefer non-sorted PSI file if both exist; both have equivalent fields, but the
# exact *.leafcutter.PSI.tsv path is the canonical source used here.
psi_files <- psi_files[!grepl("\\.leafcutter\\.PSI\\.sorted\\.tsv$", psi_files)]

sample_from_path <- function(p) {
  # .../Processed/SAMPLE/leafcutter/psi/SAMPLE.leafcutter.PSI.tsv
  parts <- strsplit(normalizePath(p, mustWork = FALSE), .Platform$file.sep, fixed = TRUE)[[1]]
  k <- which(parts == "Processed")
  if (!length(k) || k[length(k)] == length(parts)) return(NA_character_)
  parts[k[length(k)] + 1]
}

psi_index <- data.frame(
  sample_dir = vapply(psi_files, sample_from_path, character(1)),
  psi_file = psi_files,
  stringsAsFactors = FALSE
)

psi_index <- psi_index[!is.na(psi_index$sample_dir), , drop = FALSE]

meta <- merge(meta, psi_index, by = "sample_dir", all.x = TRUE)

missing_psi <- meta[is.na(meta$psi_file), c(
  "externalsampleid", "externalsubjectid", "tissue", "sample_dir", "GT", "group"
), drop = FALSE]

if (nrow(missing_psi)) {
  write.table(
    missing_psi,
    file = file.path(a$outdir, paste0(a$prefix, ".missing_psi.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

meta <- meta[!is.na(meta$psi_file), , drop = FALSE]

if (!nrow(meta)) {
  stop("None of the genotype-matched metadata samples had a LeafCutter PSI file.")
}

# Parse anchor: chr19:-:17641556-17642845
anchor_parts <- strsplit(a$anchor, ":", fixed = TRUE)[[1]]
if (length(anchor_parts) != 3) {
  stop("Anchor must look like chr19:-:17641556-17642845")
}

anchor_chr <- anchor_parts[1]
anchor_strand <- anchor_parts[2]
anchor_se <- strsplit(anchor_parts[3], "-", fixed = TRUE)[[1]]

if (length(anchor_se) != 2) {
  stop("Could not parse anchor start-end: ", anchor_parts[3])
}

anchor_start <- as.numeric(anchor_se[1])
anchor_end <- as.numeric(anchor_se[2])

if (is.na(anchor_start) || is.na(anchor_end)) {
  stop("Anchor coordinates are not numeric.")
}

read_psi <- function(path) {
  x <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      warning("Could not read ", path, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(x)) return(NULL)

  need <- c("chrom", "strand", "start", "end", "reads", "cluster_reads", "PSI")
  if (!all(need %in% colnames(x))) {
    warning("Skipping malformed PSI file: ", path)
    return(NULL)
  }

  x$start <- as.numeric(x$start)
  x$end <- as.numeric(x$end)
  x$reads <- as.numeric(x$reads)
  x
}

# -------------------------------------------------------------------------
# Event discovery
#
# Stable identity is chr:start:end:strand. We intentionally ignore sample-local
# clu_N identifiers. Starting from the sQTL anchor junction, recover the union
# across subjects of all refined LeafCutter junctions that share either the
# anchor donor/start site or anchor acceptor/end site.
#
# For UNC13A this recovers:
#   17641556-17642414
#   17641556-17642845  (anchor)
#   17642541-17642845
# -------------------------------------------------------------------------
cat("Discovering anchor-connected splice junctions...\n")

event_rows <- list()

for (i in seq_len(nrow(meta))) {
  x <- read_psi(meta$psi_file[i])
  if (is.null(x)) next

  hit <- x[
    x$chrom == anchor_chr &
      x$strand == anchor_strand &
      (x$start == anchor_start | x$end == anchor_end),
    c("chrom", "strand", "start", "end"),
    drop = FALSE
  ]

  if (nrow(hit)) event_rows[[length(event_rows) + 1]] <- hit
}

# Always retain the anchor itself.
event_rows[[length(event_rows) + 1]] <- data.frame(
  chrom = anchor_chr,
  strand = anchor_strand,
  start = anchor_start,
  end = anchor_end,
  stringsAsFactors = FALSE
)

event <- unique(do.call(rbind, event_rows))
event <- event[order(event$start, event$end), , drop = FALSE]
event$junction_id <- paste(event$chrom, event$strand,
                           paste0(event$start, "-", event$end), sep = ":")

if (!nrow(event)) {
  stop("No junctions discovered for anchor ", a$anchor)
}

cat("Discovered ", nrow(event), " anchor-connected junction(s):\n", sep = "")
cat(paste0("  ", event$junction_id, collapse = "\n"), "\n")

# -------------------------------------------------------------------------
# Build sample x junction read-count matrix.
# Missing event junctions in an otherwise valid sample are treated as 0 reads.
# Samples with total event reads == 0 are excluded from group means.
# -------------------------------------------------------------------------
count_mat <- matrix(
  0,
  nrow = nrow(meta),
  ncol = nrow(event),
  dimnames = list(meta$externalsampleid, event$junction_id)
)

for (i in seq_len(nrow(meta))) {
  x <- read_psi(meta$psi_file[i])
  if (is.null(x)) next

  x <- x[x$chrom == anchor_chr & x$strand == anchor_strand, , drop = FALSE]
  if (!nrow(x)) next

  xid <- paste(x$chrom, x$strand, paste0(x$start, "-", x$end), sep = ":")
  m <- match(event$junction_id, xid)
  ok <- !is.na(m)
  count_mat[i, ok] <- x$reads[m[ok]]
}

event_total <- rowSums(count_mat, na.rm = TRUE)
usable <- event_total > 0

sample_qc <- data.frame(
  externalsampleid = meta$externalsampleid,
  externalsubjectid = meta$externalsubjectid,
  tissue = meta$tissue,
  sample_dir = meta$sample_dir,
  GT = meta$GT,
  group = as.character(meta$group),
  event_total_reads = event_total,
  usable_event = usable,
  psi_file = meta$psi_file,
  stringsAsFactors = FALSE
)

write.table(
  sample_qc,
  file = file.path(a$outdir, paste0(a$prefix, ".samples.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (!any(usable)) {
  stop("No selected samples had reads supporting the anchor-connected event.")
}

psi_mat <- count_mat
psi_mat[,] <- NA_real_
psi_mat[usable, ] <- count_mat[usable, , drop = FALSE] / event_total[usable]

group_mean <- matrix(
  NA_real_,
  nrow = length(group_levels),
  ncol = ncol(psi_mat),
  dimnames = list(group_levels, colnames(psi_mat))
)

group_n <- setNames(integer(length(group_levels)), group_levels)

for (g in group_levels) {
  idx <- usable & as.character(meta$group) == g
  group_n[g] <- sum(idx)

  if (any(idx)) {
    group_mean[g, ] <- colMeans(psi_mat[idx, , drop = FALSE], na.rm = TRUE)
  }
}

# Diagnostic long table.
usage_long <- do.call(rbind, lapply(group_levels, function(g) {
  data.frame(
    group = g,
    n = group_n[g],
    junction_id = event$junction_id,
    chrom = event$chrom,
    strand = event$strand,
    start = event$start,
    end = event$end,
    mean_PSI = as.numeric(group_mean[g, ]),
    stringsAsFactors = FALSE
  )
}))

# -------------------------------------------------------------------------
# GTF annotation and local exon model.
# -------------------------------------------------------------------------
cat("Loading local exon annotation...\n")

region_min <- max(1, min(event$start) - 2000)
region_max <- max(event$end) + 2000

gtf_cmd <- sprintf(
  "awk -F '\\t' 'BEGIN{OFS=\"\\t\"} $1==\"%s\" && $3==\"exon\" && $4<=%d && $5>=%d {print}' %s",
  anchor_chr,
  as.integer(region_max),
  as.integer(region_min),
  shQuote(a$gtf)
)

gtf <- tryCatch(
  read.delim(
    pipe(gtf_cmd),
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  ),
  error = function(e) data.frame()
)

if (nrow(gtf)) {
  colnames(gtf)[1:9] <- c(
    "chr", "source", "feature", "start", "end", "score",
    "strand", "frame", "attributes"
  )

  extract_attr <- function(x, key) {
    pat <- paste0(key, ' "([^"]+)"')
    m <- regexec(pat, x, perl = TRUE)
    z <- regmatches(x, m)
    vapply(z, function(y) if (length(y) >= 2) y[2] else NA_character_, character(1))
  }

  gtf$gene_name <- extract_attr(gtf$attributes, "gene_name")
  gtf$gene_id <- extract_attr(gtf$attributes, "gene_id")
  gtf$transcript_id <- extract_attr(gtf$attributes, "transcript_id")

  gene_exons <- gtf[gtf$gene_name == a$gene, , drop = FALSE]
  if (!nrow(gene_exons)) {
    warning("No local GTF exons found for gene ", a$gene,
            "; using all local exons in the plotting region.")
    gene_exons <- gtf
  }
} else {
  gene_exons <- data.frame()
}

# Infer exact annotated introns from consecutive exons within transcripts.
annotated_keys <- character(0)

if (nrow(gene_exons) && any(!is.na(gene_exons$transcript_id))) {
  txs <- split(gene_exons, gene_exons$transcript_id)

  for (tx in txs) {
    tx <- tx[order(tx$start, tx$end), , drop = FALSE]
    if (nrow(tx) < 2) next

    left <- tx$end[-nrow(tx)]
    right <- tx$start[-1]
    annotated_keys <- c(
      annotated_keys,
      paste(left, right, sep = ":")
    )
  }

  annotated_keys <- unique(annotated_keys)
}

event$key <- paste(event$start, event$end, sep = ":")
event$verdict <- ifelse(event$key %in% annotated_keys, "annotated", "cryptic")
event$letter <- letters[seq_len(nrow(event))]

usage_long$verdict <- event$verdict[match(usage_long$junction_id, event$junction_id)]
usage_long$letter <- event$letter[match(usage_long$junction_id, event$junction_id)]

write.table(
  usage_long,
  file = file.path(a$outdir, paste0(a$prefix, ".junctions.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -------------------------------------------------------------------------
# LeafViz-style coordinate compression.
# -------------------------------------------------------------------------
s <- sort(unique(c(event$start, event$end)))
if (length(s) < 2) stop("Not enough unique splice sites to plot.")

length_transform <- function(g) log(g + 1)
d <- diff(s)
coords <- c(0, cumsum(length_transform(d)))
names(coords) <- as.character(s)

invert_mapping <- function(pos) {
  if (pos <= min(s)) return(coords[1])
  if (pos >= max(s)) return(coords[length(coords)])
  if (as.character(pos) %in% names(coords)) return(coords[as.character(pos)])

  w <- which(pos > s[-length(s)] & pos < s[-1])
  if (length(w) != 1) return(NA_real_)

  coords[w] +
    (coords[w + 1] - coords[w]) *
    (pos - s[w]) / (s[w + 1] - s[w])
}

total_length <- max(coords)
xlim_vals <- c(-0.05 * total_length, 1.05 * total_length)

# Exons touching event splice sites, mimicking LeafViz's culling logic.
exon_df <- data.frame()

if (nrow(gene_exons)) {
  ex <- unique(gene_exons[, c("start", "end", "gene_name", "strand"), drop = FALSE])

  ex <- ex[
    (ex$end %in% event$start | ex$start %in% event$end) &
      (
        (ex$end - ex$start) <= 500 |
          ex$end == min(event$start) |
          ex$start == max(event$end)
      ),
    ,
    drop = FALSE
  ]

  if (nrow(ex)) {
    exon_df <- data.frame(
      x = vapply(ex$start, invert_mapping, numeric(1)),
      xend = vapply(ex$end, invert_mapping, numeric(1)),
      gene_name = ex$gene_name,
      strand = ex$strand,
      stringsAsFactors = FALSE
    )

    # Ensure tiny exons remain visible.
    min_exon_length <- 0.5
    too_short <- (exon_df$xend - exon_df$x) < min_exon_length
    exon_df$xend[too_short] <- exon_df$x[too_short] + min_exon_length
    exon_df <- exon_df[!duplicated(paste(exon_df$x, exon_df$xend)), , drop = FALSE]
  }
}

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------
make_panel <- function(g, show_title = FALSE, show_legend = FALSE) {
  vals <- group_mean[g, ]
  vals[is.na(vals)] <- 0

  ed <- event
  ed$mean_PSI <- as.numeric(vals)
  ed$x <- coords[as.character(ed$start)]
  ed$xend <- coords[as.character(ed$end)]
  ed$xmid <- (ed$x + ed$xend) / 2
  ed$span <- pmax(ed$xend - ed$x, 0.01)

  # Alternate arcs above/below exactly as the original LeafViz plot does.
  ed$side <- ifelse(seq_len(nrow(ed)) %% 2 == 1, 1, -1)
  ed$ymid <- ed$side * (ed$span^0.65 / 2 + 0.25)

  # Original LeafViz effectively scales thickness ~ PSI^2.
  ed$line_size <- 0.25 + 9.75 * (ed$mean_PSI^2)
  ed$label <- paste0(
    format(ed$mean_PSI, digits = 2, nsmall = 2, scientific = FALSE),
    "  ", ed$letter
  )

  ymax <- max(abs(ed$ymid), na.rm = TRUE) * 1.35
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  p <- ggplot() +
    theme_bw(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      panel.border = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.margin = margin(5, 10, 5, 10)
    ) +
    geom_hline(yintercept = 0, colour = "white", size = 3) +
    geom_hline(yintercept = 0, colour = "black", size = 0.6)

  # Draw each arc as two curves meeting at the midpoint, following LeafViz.
  for (i in seq_len(nrow(ed))) {
    curvature <- if (ed$side[i] > 0) -0.1 else 0.1

    seg1 <- data.frame(
      x = ed$x[i], xend = ed$xmid[i],
      y = 0, yend = ed$ymid[i],
      verdict = ed$verdict[i]
    )
    seg2 <- data.frame(
      x = ed$xmid[i], xend = ed$xend[i],
      y = ed$ymid[i], yend = 0,
      verdict = ed$verdict[i]
    )

    p <- p +
      geom_curve(
        data = seg1,
        aes(x = x, xend = xend, y = y, yend = yend, colour = verdict),
        curvature = curvature,
        angle = 90,
        lineend = "round",
        size = ed$line_size[i],
        inherit.aes = FALSE
      ) +
      geom_curve(
        data = seg2,
        aes(x = x, xend = xend, y = y, yend = yend, colour = verdict),
        curvature = curvature,
        angle = 90,
        lineend = "round",
        size = ed$line_size[i],
        inherit.aes = FALSE
      ) +
      geom_label(
        data = ed[i, , drop = FALSE],
        aes(x = xmid, y = 0.92 * ymid, label = label),
        size = 3.2,
        label.size = NA,
        fill = "white",
        colour = "black",
        inherit.aes = FALSE
      )
  }

  if (nrow(exon_df)) {
    p <- p +
      geom_segment(
        data = exon_df,
        aes(x = x, xend = xend, y = 0, yend = 0),
        size = 6,
        colour = "black",
        inherit.aes = FALSE
      )

    gene_strand <- unique(exon_df$strand)
    gene_strand <- gene_strand[!is.na(gene_strand)]
    if (length(gene_strand) == 1) {
      # Add simple strand arrows in long gaps between visible exons.
      exs <- exon_df[order(exon_df$x), , drop = FALSE]
      if (nrow(exs) >= 2) {
        for (j in seq_len(nrow(exs) - 1)) {
          gap_start <- exs$xend[j]
          gap_end <- exs$x[j + 1]
          if ((gap_end - gap_start) > 0.025 * total_length) {
            mid <- (gap_start + gap_end) / 2
            if (gene_strand == "+") {
              xa <- gap_start
              xb <- mid
            } else {
              xa <- mid
              xb <- gap_end
            }

            p <- p + geom_segment(
              aes(x = xa, xend = xb, y = 0, yend = 0),
              arrow = arrow(
                ends = if (gene_strand == "+") "last" else "first",
                type = "open",
                angle = 30,
                length = unit(0.08, "inches")
              ),
              size = 0.7,
              colour = "black",
              inherit.aes = FALSE
            )
          }
        }
      }
    }
  }

  # LeafViz colors: annotated red, unannotated/cryptic pink.
  p <- p +
    scale_colour_manual(
      values = c(annotated = "red", cryptic = "pink"),
      breaks = c("annotated", "cryptic"),
      labels = c("annotated", "cryptic"),
      name = NULL,
      drop = FALSE
    ) +
    coord_cartesian(xlim = xlim_vals, ylim = c(-ymax, ymax), clip = "off") +
    ylab(sprintf("%s (n=%d)", g, group_n[g]))

  if (show_title) {
    variant_label <- paste0(
      a[["variant-chr"]], ":", a[["variant-pos"]],
      " ", a[["variant-ref"]], ">", a[["variant-alt"]]
    )
    if (!is.na(a[["variant-id"]]) && nzchar(a[["variant-id"]]) && a[["variant-id"]] != ".") {
      variant_label <- paste0(a[["variant-id"]], " | ", variant_label)
    }

    p <- p +
      ggtitle(
        a$gene,
        subtitle = paste0(a$anchor, " | ", variant_label)
      )
  }

  p
}

plots <- list()
for (i in seq_along(group_levels)) {
  g <- group_levels[i]
  plots[[i]] <- make_panel(
    g,
    show_title = (i == 1),
    show_legend = (i == length(group_levels))
  )
}

# Junction legend/table grob.
table_df <- event[, c("letter", "chrom", "strand", "start", "end", "verdict"), drop = FALSE]
colnames(table_df) <- c("ID", "chr", "strand", "start", "end", "annotation")

table_grob <- gridExtra::tableGrob(
  table_df,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 9,
    core = list(bg_params = list(fill = c("white", "grey95"), col = NA)),
    colhead = list(fg_params = list(fontface = "bold"))
  )
)

combined <- gridExtra::arrangeGrob(
  grobs = c(plots, list(table_grob)),
  ncol = 1,
  heights = c(rep(1, length(plots)), 0.8)
)

pdf_file <- file.path(a$outdir, paste0(a$prefix, ".pdf"))
png_file <- file.path(a$outdir, paste0(a$prefix, ".png"))
svg_file <- file.path(a$outdir, paste0(a$prefix, ".svg"))

ggsave(pdf_file, combined, width = 11, height = 10, units = "in")
ggsave(png_file, combined, width = 11, height = 10, units = "in", dpi = 300)
ggsave(svg_file, combined, width = 11, height = 10, units = "in")

cat("\nGenotype sample counts with event coverage:\n")
print(group_n)

cat("\nOutputs:\n")
cat("  ", pdf_file, "\n", sep = "")
cat("  ", png_file, "\n", sep = "")
cat("  ", svg_file, "\n", sep = "")
cat("  ", file.path(a$outdir, paste0(a$prefix, ".samples.tsv")), "\n", sep = "")
cat("  ", file.path(a$outdir, paste0(a$prefix, ".junctions.tsv")), "\n", sep = "")
