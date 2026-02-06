#!/usr/bin/env Rscript
# _bam2annotation.r
# prettier donut plot with external labels + legend + sample name in title

if(!require('ggplot2')) install.packages('ggplot2', repos='https://cloud.r-project.org')
if(!require('scales')) install.packages('scales', repos='https://cloud.r-project.org')
library(ggplot2)
library(scales)

args <- commandArgs(TRUE)
if(length(args) < 2) {
  stop("Usage: Rscript bam2annotation_pretty.r <Aligned.sortedByCoord.out.summary.txt> <output.pdf>")
}

input_file <- args[1]
pdf_file <- args[2]

if(!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

# infer sample name from directory
sample_name <- basename(dirname(normalizePath(input_file)))

# read summary
df <- read.table(input_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# extract values
exons <- df$exons[1]
introns <- df$introns[1]
rRNA <- df$rRNA[1]
mtRNA <- df$mtRNA[1]
intergenic_near_genes <- df$intergenic_near_genes[1]
intergenic_not_near_genes <- df$intergenic_not_near_genes[1]

# build tidy plotting table
plot_df <- data.frame(
  category=c(
    "Exons",
    "Introns",
    "Intergenic (near genes)",
    "Intergenic (not near genes)",
    "mtRNA",
    "rRNA"
  ),
  reads=c(
    exons,
    introns,
    intergenic_near_genes,
    intergenic_not_near_genes,
    mtRNA,
    rRNA
  )
)

plot_df$frac <- plot_df$reads / sum(plot_df$reads)
plot_df$label <- paste0(plot_df$category, "  ", percent(plot_df$frac, accuracy=0.1))

# color palette (colorblind-safe)
cols <- c(
  "Exons"="#4E79A7",
  "Introns"="#E15759",
  "Intergenic (near genes)"="#F28E2B",
  "Intergenic (not near genes)"="#EDC948",
  "mtRNA"="#59A14F",
  "rRNA"="#76B7B2"
)

pdf(pdf_file, width=9, height=7)

ggplot(plot_df, aes(x=2, y=reads, fill=category)) +
  geom_col(width=1, color="white") +
  coord_polar(theta="y") +
  xlim(0.5, 2.5) +
  scale_fill_manual(values=cols) +
  geom_text(
    aes(x=2, label=label),
    position=position_stack(vjust=0.5),
    hjust=0.5,
    size=3
  ) +
  ggtitle("Read Annotation Composition", subtitle=sample_name) +
  theme_void(base_size=12) +
  theme(
    legend.position="right",
    legend.title=element_blank(),
    plot.title=element_text(face="bold", hjust=0.5),
    plot.subtitle=element_text(face="italic", hjust=0.5)
  )

dev.off()
