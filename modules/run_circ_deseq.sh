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
library(ggrepel) # For better label placement

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

# --- Step 6: Volcano Plot ---
plot_data <- res %>%
  filter(!is.na(padj)) %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 0  ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0  ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Identify top 5 up and down
top_labels <- plot_data %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_max(order_by = abs(log2FoldChange), n = 5)

volcano_p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not Significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(data = top_labels, aes(label = circRNA_ID), size = 3, show.legend = FALSE) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title="Differential circRNA Expression (ALS vs Control)", y = "-log10(padj)")

ggsave(file.path(output_dir, "circRNA_volcano.png"), volcano_p, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "circRNA_volcano.svg"), volcano_p, width = 8, height = 6)
EOF
