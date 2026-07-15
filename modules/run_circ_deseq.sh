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

output_dir <- "/home/zw529/donglab/data/target_ALS"

message("Step 1: Loading and robustly cleaning metadata...")
# Use read_csv to handle embedded newlines in Cause of Death fields
meta_raw <- read_csv("~/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv", 
                     show_col_types = FALSE)

# 1. Sanitize IDs (hyphen to underscore)
# 2. Filter for valid CGND IDs to drop metadata artifacts
# 3. Keep only the first occurrence of any duplicate ID
meta_clean <- meta_raw %>%
  mutate(externalsampleid = gsub("-", "_", externalsampleid)) %>%
  filter(grepl("CGND_HRA", externalsampleid)) %>%
  distinct(externalsampleid, .keep_all = TRUE) %>%
  mutate(
    subject_group = trimws(gsub("[\r\n\t]+", " ", subject_group)),
    condition = case_when(
      grepl("ALS", subject_group) ~ "ALS",
      grepl("Control", subject_group) ~ "Control",
      TRUE ~ NA_character_
    ),
    sex = factor(sex_genotype),
    age = as.numeric(age_at_death),
    pmi = as.numeric(post_mortem_interval_in_hours)
  ) %>%
  filter(!is.na(condition), !is.na(sex), !is.na(age), !is.na(pmi))

rownames(meta_clean) <- meta_clean$externalsampleid

# Step 2: Tissue Mapping (remains same)
map_tissue <- function(t) {
  case_when(
    str_detect(t, "(?i)Motor Cortex|BA4|Lateral_motor_cortex") ~ "Motor_Cortex",
    str_detect(t, "(?i)Cervical") ~ "Cervical_Spinal_Cord",
    str_detect(t, "(?i)Lumbar|Lumbosacral") ~ "Lumbar_Spinal_Cord",
    str_detect(t, "(?i)Thoracic") ~ "Thoracic_Spinal_Cord",
    str_detect(t, "(?i)Frontal") ~ "Frontal_Cortex",
    str_detect(t, "(?i)Cerebellum") ~ "Cerebellum",
    str_detect(t, "(?i)Temporal") ~ "Temporal_Cortex",
    str_detect(t, "(?i)Medulla") ~ "Medulla",
    str_detect(t, "(?i)Pons") ~ "Pons",
    str_detect(t, "(?i)Hippocampus") ~ "Hippocampus",
    TRUE ~ "Other_Unspecified"
  )
}

meta_clean <- meta_clean %>%
  mutate(mapped_tissue = map_tissue(tissue)) %>%
  filter(mapped_tissue != "Other_Unspecified")

# Step 3: Crawl (updated for the sanitized underscores)
message("Step 3: Parsing files...")
circ_list <- list()
for (i in 1:nrow(meta_clean)) {
  s_id <- meta_clean$externalsampleid[i]
  f_path <- file.path("~/donglab/data/target_ALS", meta_clean$mapped_tissue[i], "RNAseq/Processed", s_id, "circularRNA_known_circ_percentage.txt")
  
  if (file.exists(f_path)) {
    data <- read.delim(f_path, header = TRUE)
    circ_list[[s_id]] <- data %>%
      mutate(circ_id = paste0(chrom, ":", start, "-", end, ":", strand)) %>%
      group_by(circ_id) %>% summarise(readNumber = sum(readNumber), .groups = 'drop')
  }
}

# Merge and run DESeq2
count_matrix <- purrr::reduce(circ_list, full_join, by = "circ_id") %>% replace(is.na(.), 0) %>% column_to_rownames("circ_id")
common <- intersect(colnames(count_matrix), rownames(meta_clean))
dds <- DESeqDataSetFromMatrix(as.matrix(count_matrix[, common]), meta_clean[common, ], ~ sex + age + pmi + mapped_tissue + condition)

dds <- DESeq(dds[rowSums(counts(dds) >= 5) >= 10,])
res <- as.data.frame(results(dds, contrast = c("condition", "ALS", "Control"))) %>% rownames_to_column("circRNA_ID")

write.csv(res, file.path(output_dir, "DE_circRNAs.csv"), row.names = FALSE)

# Volcano Plot
plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = pvalue < 0.05)) +
  geom_point(alpha = 0.5) + theme_minimal() + labs(title="ALS vs Control circRNA")
ggsave(file.path(output_dir, "volcano.png"), plot)
EOF
