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

# --- Step 1: Metadata ---
meta_raw <- read_csv("~/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv", show_col_types = FALSE)
meta_clean <- meta_raw %>%
  mutate(externalsampleid = gsub("-", "_", externalsampleid)) %>%
  filter(grepl("CGND_HRA", externalsampleid)) %>%
  distinct(externalsampleid, .keep_all = TRUE) %>%
  mutate(condition = if_else(grepl("ALS", subject_group), "ALS", "Control")) %>%
  filter(!is.na(age_at_death), !is.na(post_mortem_interval_in_hours))

# --- Step 2: Mapping ---
map_tissue <- function(t) {
  case_when(str_detect(t, "(?i)Motor Cortex|BA4") ~ "Motor_Cortex",
            str_detect(t, "(?i)Cervical") ~ "Cervical_Spinal_Cord",
            str_detect(t, "(?i)Lumbar") ~ "Lumbar_Spinal_Cord",
            str_detect(t, "(?i)Frontal") ~ "Frontal_Cortex",
            str_detect(t, "(?i)Cerebellum") ~ "Cerebellum",
            TRUE ~ "Other")
}
meta_clean <- meta_clean %>% mutate(mapped_tissue = map_tissue(tissue)) %>% filter(mapped_tissue != "Other")

# --- Step 3: Parsing with String Sanitization ---
circ_list <- list()
for (i in 1:nrow(meta_clean)) {
  s_id <- meta_clean$externalsampleid[i]
  f_path <- file.path("/home/zw529/donglab/data/target_ALS", meta_clean$mapped_tissue[i], "RNAseq/Processed", s_id, "circularRNA_known_circ_percentage.txt")
  
  if (file.exists(f_path)) {
    data <- read.delim(f_path, header = TRUE, stringsAsFactors = FALSE)
    if ("readNumber" %in% colnames(data)) {
      circ_list[[s_id]] <- data %>%
        # Force uniform ID formatting: remove 'chr', lowercase, trim whitespace
        mutate(circ_id = str_trim(paste0(str_remove(chrom, "chr"), ":", start, "-", end, ":", strand)),
               circ_id = str_to_lower(circ_id),
               readNumber = as.integer(readNumber)) %>%
        group_by(circ_id) %>% summarise(readNumber = sum(readNumber, na.rm = TRUE), .groups = 'drop')
    }
  }
}

# --- Step 4: Matrix Construction ---
# Use a common ID index (all circRNAs found across ANY sample)
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

dds <- DESeqDataSetFromMatrix(count_mat, colData_df, ~ sex + age + post_mortem_interval_in_hours + mapped_tissue + condition)
dds <- DESeq(dds[rowSums(counts(dds) >= 5) >= 10,])
res <- as.data.frame(results(dds, contrast = c("condition", "ALS", "Control"))) %>% rownames_to_column("circRNA_ID")

write.csv(res, file.path(output_dir, "DE_circRNAs.csv"), row.names = FALSE)
EOF
