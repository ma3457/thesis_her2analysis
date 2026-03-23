# run_gyn_HER2_scores.R
# Project breast-derived HER2 gene sets into CPTAC ovarian and endometrial proteome data
# Main goal here is to see whether the breast HER2 program shows up in gyn tumors
# and whether those projected scores track with ERBB2 protein abundance

library(tidyverse)

# --------------------------------------------------
# set working directory
# --------------------------------------------------
setwd("/Users/maya.anand/Desktop/Thesis/GYN DATA")

# --------------------------------------------------
# file paths
# --------------------------------------------------
ovarian_file <- "/Users/maya.anand/Desktop/Thesis/GYN DATA/CPTAC_OVARIAN_JHU_PNNL_FINAL_imputed_annotated (1).csv"
ucec_file    <- "/Users/maya.anand/Desktop/Thesis/GYN DATA/CPTAC_UCEC_FINAL_imputed_annotated (1).csv"

# prioritized HER2 gene sets from Aim 2
p76_file <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv"
p25_file <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_25_gene_list.csv"

# --------------------------------------------------
# read in data
# --------------------------------------------------
ovarian <- read_csv(ovarian_file, show_col_types = FALSE)
ucec    <- read_csv(ucec_file, show_col_types = FALSE)

p76_raw <- read_csv(p76_file, show_col_types = FALSE)
p25_raw <- read_csv(p25_file, show_col_types = FALSE)

# --------------------------------------------------
# grab gene sets
# your files use gene_symbol, not Gene
# --------------------------------------------------
cat("\nColumns in P76 file:\n")
print(colnames(p76_raw))

cat("\nColumns in P25 file:\n")
print(colnames(p25_raw))

if (!"gene_symbol" %in% colnames(p76_raw)) {
  stop("P76 file does not have a column named 'gene_symbol'. Check the file and column names.")
}

if (!"gene_symbol" %in% colnames(p25_raw)) {
  stop("P25 file does not have a column named 'gene_symbol'. Check the file and column names.")
}

P76 <- p76_raw %>%
  pull(gene_symbol) %>%
  na.omit() %>%
  str_trim() %>%
  unique()

P25 <- p25_raw %>%
  pull(gene_symbol) %>%
  na.omit() %>%
  str_trim() %>%
  unique()

cat("\nP76 genes loaded:", length(P76), "\n")
cat("P25 genes loaded:", length(P25), "\n")

# --------------------------------------------------
# keep only the columns we actually need
# ovarian uses TCGA sample names
# ucec has tumor and normal, so keeping tumor only
# --------------------------------------------------
ovarian_clean <- ovarian %>%
  select(Gene, starts_with("TCGA"))

ucec_clean <- ucec %>%
  select(Gene, contains("_Tumor"))

# --------------------------------------------------
# clean up gene column
# just making sure names are trimmed and not empty
# --------------------------------------------------
ovarian_clean <- ovarian_clean %>%
  mutate(Gene = str_trim(Gene))

ucec_clean <- ucec_clean %>%
  mutate(Gene = str_trim(Gene))

# --------------------------------------------------
# collapse duplicate genes
# proteomics files can have repeated gene symbols
# averaging them so we end up with one row per gene
# --------------------------------------------------
ovarian_clean <- ovarian_clean %>%
  filter(!is.na(Gene), Gene != "") %>%
  group_by(Gene) %>%
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

ucec_clean <- ucec_clean %>%
  filter(!is.na(Gene), Gene != "") %>%
  group_by(Gene) %>%
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# --------------------------------------------------
# basic checks
# --------------------------------------------------
cat("\nChecking for ERBB2...\n")
cat("ERBB2 in ovarian:", "ERBB2" %in% ovarian_clean$Gene, "\n")
cat("ERBB2 in UCEC:", "ERBB2" %in% ucec_clean$Gene, "\n")

cat("\nDimensions after cleanup...\n")
cat("Ovarian:", nrow(ovarian_clean), "genes x", ncol(ovarian_clean) - 1, "samples\n")
cat("UCEC:", nrow(ucec_clean), "genes x", ncol(ucec_clean) - 1, "tumor samples\n")

# --------------------------------------------------
# check how many signature genes are detectable
# this matters because proteome coverage is always lower than RNA
# --------------------------------------------------
cat("\nGene set overlap counts...\n")
cat("P76 in ovarian:", sum(P76 %in% ovarian_clean$Gene), "\n")
cat("P76 in UCEC:", sum(P76 %in% ucec_clean$Gene), "\n")
cat("P25 in ovarian:", sum(P25 %in% ovarian_clean$Gene), "\n")
cat("P25 in UCEC:", sum(P25 %in% ucec_clean$Gene), "\n")

# --------------------------------------------------
# subset to only detected genes from each signature
# --------------------------------------------------
ovarian_P76 <- ovarian_clean %>% filter(Gene %in% P76)
ovarian_P25 <- ovarian_clean %>% filter(Gene %in% P25)

ucec_P76 <- ucec_clean %>% filter(Gene %in% P76)
ucec_P25 <- ucec_clean %>% filter(Gene %in% P25)

# --------------------------------------------------
# helper function:
# z-score each gene across samples, then average those z-scores by sample
# this gives a per-sample HER2-related signature score
# --------------------------------------------------
compute_signature_score <- function(df) {
  if (nrow(df) == 0) {
    stop("No genes found in dataset for this signature.")
  }
  
  mat <- df %>%
    tibble::column_to_rownames("Gene") %>%
    as.matrix()
  
  # z-score each gene across samples
  mat_z <- t(scale(t(mat)))
  
  # if a gene has zero variance, scale gives NA
  # setting those to 0 keeps them from breaking the score
  mat_z[is.na(mat_z)] <- 0
  
  # average z-scored signal across genes for each sample
  scores <- colMeans(mat_z, na.rm = TRUE)
  return(scores)
}

# --------------------------------------------------
# calculate scores
# --------------------------------------------------
ovarian_P76_scores <- compute_signature_score(ovarian_P76)
ovarian_P25_scores <- compute_signature_score(ovarian_P25)

ucec_P76_scores <- compute_signature_score(ucec_P76)
ucec_P25_scores <- compute_signature_score(ucec_P25)

# --------------------------------------------------
# get ERBB2 protein abundance
# --------------------------------------------------
get_erbb2 <- function(df) {
  erbb2_row <- df %>% filter(Gene == "ERBB2")
  
  if (nrow(erbb2_row) == 0) {
    stop("ERBB2 not found in dataset.")
  }
  
  erbb2_vals <- erbb2_row %>%
    select(-Gene) %>%
    as.numeric()
  
  names(erbb2_vals) <- colnames(erbb2_row %>% select(-Gene))
  return(erbb2_vals)
}

ovarian_ERBB2 <- get_erbb2(ovarian_clean)
ucec_ERBB2 <- get_erbb2(ucec_clean)

# --------------------------------------------------
# build output data frames
# using sample names to line things up correctly
# --------------------------------------------------
ovarian_plot_df <- tibble(
  Sample = names(ovarian_P76_scores),
  P76_score = as.numeric(ovarian_P76_scores),
  P25_score = as.numeric(ovarian_P25_scores[names(ovarian_P76_scores)]),
  ERBB2 = as.numeric(ovarian_ERBB2[names(ovarian_P76_scores)])
)

ucec_plot_df <- tibble(
  Sample = names(ucec_P76_scores),
  P76_score = as.numeric(ucec_P76_scores),
  P25_score = as.numeric(ucec_P25_scores[names(ucec_P76_scores)]),
  ERBB2 = as.numeric(ucec_ERBB2[names(ucec_P76_scores)])
)

# --------------------------------------------------
# save score tables
# --------------------------------------------------
write_csv(ovarian_plot_df, "ovarian_HER2_scores.csv")
write_csv(ucec_plot_df, "ucec_HER2_scores.csv")

# --------------------------------------------------
# save detected genes too
# useful for methods/results wording later
# --------------------------------------------------
write_csv(tibble(Gene = ovarian_P76$Gene), "ovarian_detected_P76_genes.csv")
write_csv(tibble(Gene = ovarian_P25$Gene), "ovarian_detected_P25_genes.csv")
write_csv(tibble(Gene = ucec_P76$Gene), "ucec_detected_P76_genes.csv")
write_csv(tibble(Gene = ucec_P25$Gene), "ucec_detected_P25_genes.csv")

# --------------------------------------------------
# correlation checks
# first pass look at whether ERBB2 tracks with the projected signature
# --------------------------------------------------
ovarian_cor_p76 <- cor.test(ovarian_plot_df$ERBB2, ovarian_plot_df$P76_score, method = "spearman")
ovarian_cor_p25 <- cor.test(ovarian_plot_df$ERBB2, ovarian_plot_df$P25_score, method = "spearman")

ucec_cor_p76 <- cor.test(ucec_plot_df$ERBB2, ucec_plot_df$P76_score, method = "spearman")
ucec_cor_p25 <- cor.test(ucec_plot_df$ERBB2, ucec_plot_df$P25_score, method = "spearman")

cat("\nSpearman correlations with ERBB2...\n")
cat("Ovarian P76 rho:", round(ovarian_cor_p76$estimate, 3), " p =", signif(ovarian_cor_p76$p.value, 3), "\n")
cat("Ovarian P25 rho:", round(ovarian_cor_p25$estimate, 3), " p =", signif(ovarian_cor_p25$p.value, 3), "\n")
cat("UCEC P76 rho:", round(ucec_cor_p76$estimate, 3), " p =", signif(ucec_cor_p76$p.value, 3), "\n")
cat("UCEC P25 rho:", round(ucec_cor_p25$estimate, 3), " p =", signif(ucec_cor_p25$p.value, 3), "\n")

# --------------------------------------------------
# scatterplots
# keeping these simple for now as first-pass figures
# --------------------------------------------------
p1 <- ggplot(ovarian_plot_df, aes(x = ERBB2, y = P76_score)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Ovarian CPTAC: P76 HER2 score vs ERBB2 protein",
    x = "ERBB2 protein abundance",
    y = "P76 HER2 score"
  ) +
  theme_classic()

p2 <- ggplot(ovarian_plot_df, aes(x = ERBB2, y = P25_score)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Ovarian CPTAC: P25 HER2 score vs ERBB2 protein",
    x = "ERBB2 protein abundance",
    y = "P25 HER2 score"
  ) +
  theme_classic()

p3 <- ggplot(ucec_plot_df, aes(x = ERBB2, y = P76_score)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "UCEC CPTAC: P76 HER2 score vs ERBB2 protein",
    x = "ERBB2 protein abundance",
    y = "P76 HER2 score"
  ) +
  theme_classic()

p4 <- ggplot(ucec_plot_df, aes(x = ERBB2, y = P25_score)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "UCEC CPTAC: P25 HER2 score vs ERBB2 protein",
    x = "ERBB2 protein abundance",
    y = "P25 HER2 score"
  ) +
  theme_classic()

ggsave("ovarian_P76_vs_ERBB2.png", p1, width = 7, height = 5, dpi = 300)
ggsave("ovarian_P25_vs_ERBB2.png", p2, width = 7, height = 5, dpi = 300)
ggsave("ucec_P76_vs_ERBB2.png", p3, width = 7, height = 5, dpi = 300)
ggsave("ucec_P25_vs_ERBB2.png", p4, width = 7, height = 5, dpi = 300)

# --------------------------------------------------
# top scoring tumors
# these are probably the ones worth checking first
# --------------------------------------------------
ovarian_top <- ovarian_plot_df %>%
  arrange(desc(P76_score)) %>%
  slice(1:10)

ucec_top <- ucec_plot_df %>%
  arrange(desc(P76_score)) %>%
  slice(1:10)

write_csv(ovarian_top, "ovarian_top_P76_samples.csv")
write_csv(ucec_top, "ucec_top_P76_samples.csv")

cat("\nTop ovarian samples by P76 score saved.\n")
cat("Top UCEC samples by P76 score saved.\n")

cat("\nDone.\n")