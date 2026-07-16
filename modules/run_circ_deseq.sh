#!/bin/bash
#SBATCH --job-name=circ_deseq
#SBATCH --output=/home/zw529/donglab/data/target_ALS/circ_deseq_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/circ_deseq_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

module load R

mkdir -p ~/donglab/data/target_ALS

Rscript - <<'EOF'
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(readr)
library(ggrepel)

output_dir <- "/home/zw529/donglab/data/target_ALS"

# --- Step 1 & 2: Metadata and Mapping ---
meta_raw <- read_csv("~/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv", show_col_types = FALSE)
meta_clean <- meta_raw %>%
  mutate(externalsampleid = gsub("-", "_", externalsampleid)) %>%
  filter(grepl("CGND_HRA", externalsampleid)) %>%
  distinct(externalsampleid, .keep_all = TRUE) %>%
  mutate(condition = if_else(grepl("ALS", subject_group), "ALS", "Control")) %>%
  rename(age = age_at_death, pmi = post_mortem_interval_in_hours) %>%
  filter(!is.na(age), !is.na(pmi))

map_tissue <- function(t) {
  case_when(str_detect(t, "(?i)Motor Cortex|BA4") ~ "Motor_Cortex",
            str_detect(t, "(?i)Cervical") ~ "Cervical_Spinal_Cord",
            str_detect(t, "(?i)Lumbar") ~ "Lumbar_Spinal_Cord",
            str_detect(t, "(?i)Frontal") ~ "Frontal_Cortex",
            str_detect(t, "(?i)Cerebellum") ~ "Cerebellum",
            TRUE ~ "Other")
}
meta_clean <- meta_clean %>% mutate(mapped_tissue = map_tissue(tissue)) %>% filter(mapped_tissue != "Other")

# --- Step 3 & 4: Parsing and Matrix ---
circ_list <- list()
for (i in 1:nrow(meta_clean)) {
  s_id <- meta_clean$externalsampleid[i]
  f_path <- file.path("/home/zw529/donglab/data/target_ALS", meta_clean$mapped_tissue[i], "RNAseq/Processed", s_id, "circularRNA_known_circ_percentage.txt")
  if (file.exists(f_path)) {
    data <- read.delim(f_path, header = TRUE, stringsAsFactors = FALSE)
    if ("readNumber" %in% colnames(data)) {
      circ_list[[s_id]] <- data %>%
        mutate(circ_id = str_to_lower(str_trim(paste0(str_remove(chrom, "chr"), ":", start, "-", end, ":", strand))),
               readNumber = as.integer(readNumber)) %>%
        group_by(circ_id) %>% summarise(readNumber = sum(readNumber, na.rm = TRUE), .groups = 'drop')
    }
  }
}

all_ids <- unique(unlist(lapply(circ_list, function(x) x$circ_id)))
count_mat <- matrix(0, nrow = length(all_ids), ncol = length(circ_list), dimnames = list(all_ids, names(circ_list)))
for (s_id in names(circ_list)) {
  df <- circ_list[[s_id]]
  count_mat[df$circ_id, s_id] <- df$readNumber
}

# --- Step 5: DESeq2 ---
common_samples <- intersect(colnames(count_mat), meta_clean$externalsampleid)
count_mat <- count_mat[, common_samples]
colData_df <- as.data.frame(meta_clean %>% filter(externalsampleid %in% common_samples))
rownames(colData_df) <- colData_df$externalsampleid

dds <- DESeqDataSetFromMatrix(count_mat, colData_df, ~ sex + age + pmi + mapped_tissue + condition)
dds <- DESeq(dds[rowSums(counts(dds) >= 5) >= 10,])
res <- as.data.frame(results(dds, contrast = c("condition", "ALS", "Control"))) %>% rownames_to_column("circRNA_ID")

write.csv(res, file.path(output_dir, "DE_circRNAs.csv"), row.names = FALSE)

# --- Step 6: Annotation & Volcano Plot ---
# Load and parse BED for gene names
bed_file <- "~/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6"
gene_anno <- read.table(bed_file, sep = "\t", header = FALSE) %>%
  select(chrom = V1, start = V2, end = V3, info = V4) %>%
  mutate(gene_name = str_split_fixed(info, "___", 3)[,3]) %>%
  select(chrom, start, end, gene_name)

res_annotated <- res %>%
  mutate(
    temp_chr = str_remove(circRNA_ID, ":.*"),
    start = as.numeric(str_extract(str_remove(circRNA_ID, ".*:"), "\\d+")),
    chrom = paste0("chr", str_remove(temp_chr, "chr"))
  ) %>%
  left_join(gene_anno, by = c("chrom", "start")) %>%
  mutate(label = ifelse(!is.na(gene_name), paste0(gene_name, "\n(chr", circRNA_ID, ")"), paste0("chr", circRNA_ID)))

plot_data <- res_annotated %>%
  filter(!is.na(padj)) %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 0  ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0  ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

top_labels <- plot_data %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_max(order_by = abs(log2FoldChange), n = 4)

volcano_p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not Significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(data = top_labels, aes(label = label), size = 2.5, box.padding = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title="Differential circRNA Expression (ALS vs Control)", 
       color = "ALS-associated expression",
       y = "-log10(padj)")

ggsave(file.path(output_dir, "circRNA_volcano.png"), volcano_p, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "circRNA_volcano.svg"), volcano_p, width = 8, height = 6)
EOF
