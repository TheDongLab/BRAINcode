#!/bin/bash
#SBATCH --job-name=circ_deseq
#SBATCH --output=/home/zw529/donglab/data/target_ALS/circ_deseq_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/circ_deseq_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

module load R

# Create the output directory if it doesn't exist
mkdir -p ~/donglab/data/target_ALS

# Execute R using a heredoc
Rscript - <<'EOF'
# ==========================================
# R SCRIPT WRAPPED INSIDE BATCH
# ==========================================
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

output_dir <- "~/donglab/data/target_ALS"

message("Step 1: Loading and cleaning metadata...")
meta_raw <- read.csv("~/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv", 
                     stringsAsFactors = FALSE, 
                     na.strings = c("", "NA", "N/A", "unknown", "Unknown"))

# Apply cleanups and drop missing covariate rows
meta_clean <- meta_raw %>%
  mutate(
    subject_group = trimws(gsub("[\r\n\t]+", " ", subject_group)),
    subject_group = case_when(
      subject_group == "Non Neurological Control" ~ "Non-Neurological Control",
      subject_group == "ALS Spectrum MND, Other Neurological Diseases" ~ "ALS Spectrum MND, Other Neurological Disorders",
      TRUE ~ subject_group
    )
  ) %>%
  mutate(
    condition = case_match(
      subject_group,
      "ALS Spectrum MND" ~ "ALS",
      "Non-Neurological Control" ~ "Control",
      .default = NA_character_
    )
  ) %>%
  mutate(
    sex = factor(sex_genotype),
    age = as.numeric(age_at_death),
    pmi = as.numeric(post_mortem_interval_in_hours)
  ) %>%
  filter(!is.na(condition), !is.na(sex), !is.na(age), !is.na(pmi))

# ==========================================
# Step 2: Standardize Tissue Mappings
# ==========================================
map_tissue <- function(t) {
  case_when(
    str_detect(t, "(?i)Motor Cortex|BA4|Lateral_motor_cortex") ~ "Motor_Cortex",
    str_detect(t, "(?i)Cervical")                              ~ "Cervical_Spinal_Cord",
    str_detect(t, "(?i)Lumbar|Lumbosacral")                    ~ "Lumbar_Spinal_Cord",
    str_detect(t, "(?i)Thoracic")                              ~ "Thoracic_Spinal_Cord",
    str_detect(t, "(?i)Frontal")                               ~ "Frontal_Cortex",
    str_detect(t, "(?i)Cerebellum")                            ~ "Cerebellum",
    str_detect(t, "(?i)Temporal")                              ~ "Temporal_Cortex",
    str_detect(t, "(?i)Medulla")                               ~ "Medulla",
    str_detect(t, "(?i)Pons")                                  ~ "Pons",
    str_detect(t, "(?i)Hippocampus")                           ~ "Hippocampus",
    TRUE                                                       ~ "Other_Unspecified"
  )
}

meta_clean <- meta_clean %>%
  mutate(mapped_tissue = map_tissue(tissue)) %>%
  filter(mapped_tissue != "Other_Unspecified")

rownames(meta_clean) <- meta_clean$externalsampleid

# ==========================================
# Step 3: Crawl and Parse readNumber counts
# ==========================================
message("Step 3: Parsing circularRNA_known_circ_percentage.txt files...")
circ_list <- list()

for (i in 1:nrow(meta_clean)) {
  sample_id <- meta_clean$externalsampleid[i]
  sample_dir <- gsub("-", "_", sample_id) 
  tissue_dir <- meta_clean$mapped_tissue[i]
  
  file_path <- file.path("~/donglab/data/target_ALS", tissue_dir, "RNAseq/Processed", sample_dir, "circularRNA_known_circ_percentage.txt")
  
  if (file.exists(file_path)) {
    circ_data <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
    
    if (nrow(circ_data) > 0) {
      circ_data <- circ_data %>%
        mutate(circ_id = paste0(chrom, ":", start, "-", end, ":", strand)) %>%
        select(circ_id, readNumber) %>%
        group_by(circ_id) %>%
        summarise(readNumber = sum(readNumber), .groups = "drop")
      
      circ_list[[sample_id]] <- circ_data
    }
  }
}

if (length(circ_list) == 0) {
  stop("Error: No valid circularRNA_known_circ_percentage.txt files were processed.")
}

message("Merging parsed matrices...")
count_matrix <- purrr::reduce(circ_list, function(df1, df2) {
  full_join(df1, df2, by = "circ_id")
}) %>% replace(is.na(.), 0)

count_matrix <- count_matrix %>% column_to_rownames("circ_id")

common_samples <- intersect(colnames(count_matrix), rownames(meta_clean))
counts <- as.matrix(count_matrix[, common_samples])
colData <- meta_clean[common_samples, ]

# ==========================================
# Step 4: Run DESeq2
# ==========================================
message("Step 4: Setting up and running DESeq2...")
colData$condition <- factor(colData$condition, levels = c("Control", "ALS"))
colData$sex <- factor(colData$sex)
colData$mapped_tissue <- factor(colData$mapped_tissue)

# Run multi-tissue framework (adding mapped_tissue as a covariate)
design_formula <- ~ sex + age + pmi + mapped_tissue + condition

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = design_formula)

# Moderate pre-filtering
keep <- rowSums(counts(dds) >= 5) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "ALS", "Control"))

# Format results
res_df <- as.data.frame(res) %>% 
  rownames_to_column("circRNA_ID") %>% 
  arrange(pvalue)

# Save results file
output_csv <- file.path(output_dir, "DE_circRNAs_multi_tissue.csv")
write.csv(res_df, file = output_csv, row.names = FALSE)
message("Results saved to: ", output_csv)

# ==========================================
# Step 5: Generate Volcano Plot
# ==========================================
message("Step 5: Plotting Volcano Plot...")

# Determine significance groupings
plot_data <- res_df %>%
  filter(!is.na(pvalue)) %>%
  mutate(
    neg_log10_p = -log10(pvalue),
    Significance = case_when(
      pvalue < 0.05 & log2FoldChange > 0.5  ~ "Upregulated (p < 0.05)",
      pvalue < 0.05 & log2FoldChange < -0.5 ~ "Downregulated (p < 0.05)",
      TRUE                                  ~ "Not Significant"
    )
  )

volcano_p <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_p, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated (p < 0.05)" = "#d95f02", 
                                "Downregulated (p < 0.05)" = "#1f78b4", 
                                "Not Significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Differential Expression of circRNAs (ALS vs. Control)",
    subtitle = "Adjusted for Sex, Age, PMI, and Tissue Type",
    x = "log2(Fold Change)",
    y = "-log10(p-value)"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30")
  )

output_png <- file.path(output_dir, "circRNA_volcano_plot.png")
ggsave(filename = output_png, plot = volcano_p, width = 8, height = 7, dpi = 300)
message("Volcano plot successfully saved to: ", output_png)

EOF
