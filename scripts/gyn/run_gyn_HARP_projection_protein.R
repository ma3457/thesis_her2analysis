# ============================================================
# HARP candidate projection into GYN proteomes
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
})

# ------------------------------------------------------------
# FILE PATHS
# ------------------------------------------------------------
candidate_file <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/Wolf_HARP_RNA/P76_HARP_core_cotrending_candidates.csv"

ovarian_file <- "/Users/maya.anand/Desktop/Thesis/GYN DATA/CPTAC_OVARIAN_JHU_PNNL_FINAL_imputed_annotated (1).csv"
ucec_file    <- "/Users/maya.anand/Desktop/Thesis/GYN DATA/CPTAC_UCEC_FINAL_imputed_annotated (1).csv"

outdir <- "/Users/maya.anand/Desktop/Thesis/GYN DATA/HARP_projection_outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------
clean_gene <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- toupper(x)
  x[x == "" | is.na(x)] <- NA
  x
}

load_gyn_proteome <- function(file, cohort = c("ovarian", "ucec")) {
  cohort <- match.arg(cohort)
  
  df <- read_csv(file, show_col_types = FALSE)
  
  if (!"Gene" %in% colnames(df)) {
    stop(paste("Gene column not found in:", file))
  }
  
  df$Gene <- clean_gene(df$Gene)
  df <- df %>% filter(!is.na(Gene))
  
  # Sample columns:
  # Ovarian sample columns begin with TCGA.
  # UCEC file contains tumor + normal, so keep tumor only.
  if (cohort == "ovarian") {
    sample_cols <- grep("^TCGA\\.", colnames(df), value = TRUE)
  } else if (cohort == "ucec") {
    sample_cols <- grep("_Tumor_", colnames(df), value = TRUE)
  }
  
  if (length(sample_cols) == 0) {
    stop(paste("No sample columns detected for cohort:", cohort))
  }
  
  mat_df <- df[, c("Gene", sample_cols), drop = FALSE]
  
  # Convert sample columns to numeric
  for (cc in sample_cols) {
    mat_df[[cc]] <- as.numeric(mat_df[[cc]])
  }
  
  # Collapse duplicated genes by median across rows
  mat_df <- mat_df %>%
    group_by(Gene) %>%
    summarize(across(all_of(sample_cols), ~ median(.x, na.rm = TRUE)), .groups = "drop")
  
  mat <- as.data.frame(mat_df)
  rownames(mat) <- mat$Gene
  mat$Gene <- NULL
  mat <- as.matrix(mat)
  
  return(mat)
}

correlate_candidates_with_erbb2 <- function(mat, candidates_df, cohort_name) {
  # Candidate genes
  if (!"feature_id" %in% colnames(candidates_df)) {
    stop("Candidate file must contain 'feature_id'")
  }
  if (!"harp_logFC" %in% colnames(candidates_df)) {
    stop("Candidate file must contain 'harp_logFC'")
  }
  
  candidates_df <- candidates_df %>%
    mutate(feature_id = clean_gene(feature_id)) %>%
    filter(!is.na(feature_id)) %>%
    distinct(feature_id, .keep_all = TRUE)
  
  genes_all <- candidates_df$feature_id
  genes_present <- intersect(genes_all, rownames(mat))
  genes_missing <- setdiff(genes_all, rownames(mat))
  
  if (!"ERBB2" %in% rownames(mat)) {
    stop(paste("ERBB2 not found in", cohort_name, "matrix rownames"))
  }
  
  erbb2_vec <- as.numeric(mat["ERBB2", ])
  
  res <- lapply(genes_present, function(g) {
    gene_vec <- as.numeric(mat[g, ])
    ok <- is.finite(gene_vec) & is.finite(erbb2_vec)
    
    if (sum(ok) < 3) {
      return(data.frame(
        Gene = g,
        Cohort = cohort_name,
        N = sum(ok),
        Spearman_rho = NA_real_,
        Spearman_p = NA_real_
      ))
    }
    
    ct <- suppressWarnings(cor.test(
      x = gene_vec[ok],
      y = erbb2_vec[ok],
      method = "spearman",
      exact = FALSE
    ))
    
    data.frame(
      Gene = g,
      Cohort = cohort_name,
      N = sum(ok),
      Spearman_rho = unname(ct$estimate),
      Spearman_p = ct$p.value
    )
  }) %>% bind_rows()
  
  res <- res %>%
    mutate(FDR = p.adjust(Spearman_p, method = "BH")) %>%
    left_join(
      candidates_df %>%
        transmute(
          Gene = feature_id,
          Wolf_HER2_logFC = her2_logFC,
          Wolf_HARP_logFC = harp_logFC,
          Wolf_HER2_FDR = her2_FDR,
          Wolf_HARP_FDR = harp_FDR
        ),
      by = "Gene"
    ) %>%
    arrange(desc(Spearman_rho))
  
  summary_tbl <- data.frame(
    Cohort = cohort_name,
    Candidates_in_file = length(genes_all),
    Candidates_detected = length(genes_present),
    Candidates_missing = length(genes_missing),
    ERBB2_present = "ERBB2" %in% rownames(mat),
    Median_rho = median(res$Spearman_rho, na.rm = TRUE),
    Mean_rho = mean(res$Spearman_rho, na.rm = TRUE),
    Positive_rho_n = sum(res$Spearman_rho > 0, na.rm = TRUE),
    Negative_rho_n = sum(res$Spearman_rho < 0, na.rm = TRUE),
    FDR_sig_n = sum(res$FDR < 0.05, na.rm = TRUE)
  )
  
  missing_tbl <- data.frame(
    Cohort = cohort_name,
    Missing_Gene = genes_missing
  )
  
  list(results = res, summary = summary_tbl, missing = missing_tbl)
}

make_barplot <- function(res_df, cohort_name, outfile) {
  plot_df <- res_df %>%
    filter(!is.na(Spearman_rho)) %>%
    arrange(Spearman_rho) %>%
    mutate(Gene = factor(Gene, levels = Gene))
  
  p <- ggplot(plot_df, aes(x = Gene, y = Spearman_rho)) +
    geom_col() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(cohort_name, ": HArP core candidates vs ERBB2 abundance"),
      x = NULL,
      y = "Spearman rho"
    )
  
  ggsave(outfile, plot = p, width = 8, height = 6, dpi = 300)
  p
}

make_cross_context_scatter <- function(res_df, cohort_name, outfile) {
  plot_df <- res_df %>%
    filter(!is.na(Wolf_HARP_logFC), !is.na(Spearman_rho))
  
  p <- ggplot(plot_df, aes(x = Wolf_HARP_logFC, y = Spearman_rho, label = Gene)) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(cohort_name, ": Wolf HArP effect vs GYN ERBB2 association"),
      x = "Wolf HArP logFC",
      y = "Spearman rho with ERBB2 protein"
    )
  
  ggsave(outfile, plot = p, width = 6.5, height = 5.5, dpi = 300)
  p
}

# ------------------------------------------------------------
# LOAD INPUTS
# ------------------------------------------------------------
candidates_df <- read_csv(candidate_file, show_col_types = FALSE)

cat("\nCandidate file columns:\n")
print(colnames(candidates_df))

ovarian_mat <- load_gyn_proteome(ovarian_file, cohort = "ovarian")
ucec_mat    <- load_gyn_proteome(ucec_file, cohort = "ucec")

cat("\nOvarian matrix dimensions:\n")
print(dim(ovarian_mat))

cat("\nUCEC tumor-only matrix dimensions:\n")
print(dim(ucec_mat))

cat("\nERBB2 present?\n")
print(c(
  ovarian_ERBB2 = "ERBB2" %in% rownames(ovarian_mat),
  ucec_ERBB2    = "ERBB2" %in% rownames(ucec_mat)
))

# ------------------------------------------------------------
# RUN ANALYSIS
# ------------------------------------------------------------
ovarian_out <- correlate_candidates_with_erbb2(
  mat = ovarian_mat,
  candidates_df = candidates_df,
  cohort_name = "Ovarian"
)

ucec_out <- correlate_candidates_with_erbb2(
  mat = ucec_mat,
  candidates_df = candidates_df,
  cohort_name = "UCEC"
)

combined_summary <- bind_rows(ovarian_out$summary, ucec_out$summary)
combined_results <- bind_rows(ovarian_out$results, ucec_out$results)
combined_missing <- bind_rows(ovarian_out$missing, ucec_out$missing)

# ------------------------------------------------------------
# SAVE TABLES
# ------------------------------------------------------------
write_csv(ovarian_out$results, file.path(outdir, "ovarian_HArP_candidate_vs_ERBB2_results.csv"))
write_csv(ucec_out$results, file.path(outdir, "ucec_HArP_candidate_vs_ERBB2_results.csv"))
write_csv(combined_results, file.path(outdir, "combined_HArP_candidate_vs_ERBB2_results.csv"))

write_csv(ovarian_out$summary, file.path(outdir, "ovarian_HArP_candidate_summary.csv"))
write_csv(ucec_out$summary, file.path(outdir, "ucec_HArP_candidate_summary.csv"))
write_csv(combined_summary, file.path(outdir, "combined_HArP_candidate_summary.csv"))

write_csv(combined_missing, file.path(outdir, "missing_HArP_candidates_by_cohort.csv"))

# ------------------------------------------------------------
# SAVE PLOTS
# ------------------------------------------------------------
make_barplot(
  ovarian_out$results,
  "Ovarian",
  file.path(outdir, "ovarian_HArP_candidate_barplot.png")
)

make_barplot(
  ucec_out$results,
  "UCEC",
  file.path(outdir, "ucec_HArP_candidate_barplot.png")
)

make_cross_context_scatter(
  ovarian_out$results,
  "Ovarian",
  file.path(outdir, "ovarian_HArP_cross_context_scatter.png")
)

make_cross_context_scatter(
  ucec_out$results,
  "UCEC",
  file.path(outdir, "ucec_HArP_cross_context_scatter.png")
)

# ------------------------------------------------------------
# PRINT QUICK OUTPUT
# ------------------------------------------------------------
cat("\n================ SUMMARY ================\n")
print(combined_summary)

cat("\nTop Ovarian genes:\n")
print(ovarian_out$results %>% arrange(desc(Spearman_rho)) %>% head(10))

cat("\nTop UCEC genes:\n")
print(ucec_out$results %>% arrange(desc(Spearman_rho)) %>% head(10))

cat("\nOutput directory:\n")
cat(outdir, "\n")
cat("=========================================\n")
