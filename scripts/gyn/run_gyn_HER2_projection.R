# ===============================
# run_gyn_HER2_projection.R
# Project Wolf HARP-consistent HER2 candidates into
# CPTAC ovarian and UCEC proteomic cohorts
# Primary analysis: DE (ERBB2-high vs ERBB2-low)
# Secondary analysis: correlation with ERBB2 abundance
# Runs TWO Wolf-derived signatures:
#   1) P76_HARP_core_cotrending_candidates
#   2) P76_HARP_core_cotrending_candidates_UP
# ===============================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(limma)
})

# -------------------------------
# PATHS
# -------------------------------
base_dir <- "/Users/maya.anand/Desktop/Thesis/GYN DATA"

ov_fp   <- file.path(base_dir, "CPTAC_OVARIAN_JHU_PNNL_FINAL_imputed_annotated (1).csv")
ucec_fp <- file.path(base_dir, "CPTAC_UCEC_FINAL_imputed_annotated (1).csv")

sig_files <- list(
  P76_HARP_core_cotrending =
    "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/Wolf_HARP_RNA/P76_HARP_core_cotrending_candidates.csv",
  P76_HARP_core_cotrending_UP =
    "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/Wolf_HARP_RNA/P76_HARP_core_cotrending_candidates_UP.csv"
)

out_dir <- file.path(base_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(ov_fp), file.exists(ucec_fp))
stopifnot(all(file.exists(unlist(sig_files))))

# -------------------------------
# HELPERS
# -------------------------------
get_gene_col <- function(dt) {
  possible <- c("gene", "Gene", "feature_id", "feature", "gene_symbol", "GeneSymbol", "symbol")
  hit <- possible[possible %in% names(dt)][1]
  if (is.na(hit)) stop("No recognizable gene column found.")
  hit
}

load_signature <- function(fp) {
  dt <- fread(fp)
  gene_col <- get_gene_col(dt)
  unique(toupper(as.character(dt[[gene_col]])))
}

detect_gene_col <- function(dt) {
  possible <- c("gene", "Gene", "geneName", "gene_symbol", "GeneSymbol", "symbol",
                "feature_id", "feature", "Protein", "protein")
  hit <- possible[possible %in% names(dt)][1]
  if (is.na(hit)) {
    hit <- names(dt)[1]
    message("Using first column as gene column: ", hit)
  }
  hit
}

prepare_matrix <- function(fp, cohort_name) {
  dt <- fread(fp)
  
  gene_col <- detect_gene_col(dt)
  dt[[gene_col]] <- toupper(trimws(as.character(dt[[gene_col]])))
  
  dt <- dt[!is.na(get(gene_col)) & get(gene_col) != ""]
  dt <- dt[!duplicated(get(gene_col))]
  
  sample_cols <- setdiff(names(dt), gene_col)
  mat <- as.matrix(dt[, ..sample_cols])
  mode(mat) <- "numeric"
  rownames(mat) <- dt[[gene_col]]
  
  keep_rows <- rowSums(is.finite(mat)) > 0
  mat <- mat[keep_rows, , drop = FALSE]
  
  if (!"ERBB2" %in% rownames(mat)) {
    stop("ERBB2 not found in ", cohort_name, ". Check gene column / row names.")
  }
  
  mat
}

run_projection_analysis <- function(mat, cohort_name, signature_genes, signature_name, out_dir) {
  
  cat("\nRunning cohort:", cohort_name, "| signature:", signature_name, "\n")
  
  # ---------------------------
  # ERBB2 abundance and quartiles
  # ---------------------------
  erbb2 <- as.numeric(mat["ERBB2", ])
  names(erbb2) <- colnames(mat)
  
  q_low  <- unname(quantile(erbb2, 0.25, na.rm = TRUE))
  q_high <- unname(quantile(erbb2, 0.75, na.rm = TRUE))
  
  group <- ifelse(erbb2 <= q_low, "LOW",
                  ifelse(erbb2 >= q_high, "HIGH", NA))
  
  keep_samples <- !is.na(group)
  mat_q <- mat[, keep_samples, drop = FALSE]
  erbb2_q <- erbb2[keep_samples]
  group <- factor(group[keep_samples], levels = c("LOW", "HIGH"))
  
  cat("Samples used (quartiles only):", length(group), "\n")
  cat("LOW:", sum(group == "LOW"), " HIGH:", sum(group == "HIGH"), "\n")
  
  # ---------------------------
  # LIMMA DE: HIGH vs LOW
  # ---------------------------
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("LOW", "HIGH")
  
  fit <- lmFit(mat_q, design)
  cont <- makeContrasts(HIGH - LOW, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont), trend = TRUE, robust = TRUE)
  
  de_full <- as.data.table(topTable(fit2, number = Inf, sort.by = "P"))
  de_full[, gene := rownames(topTable(fit2, number = Inf, sort.by = "P"))]
  
  # ---------------------------
  # Restrict to signature genes
  # ---------------------------
  de_sig <- de_full[toupper(gene) %in% signature_genes]
  
  # ---------------------------
  # Correlation with ERBB2 across all samples
  # ---------------------------
  cor_dt <- data.table(
    gene = rownames(mat),
    rho = apply(mat, 1, function(x) {
      suppressWarnings(cor(as.numeric(x), erbb2, method = "spearman", use = "complete.obs"))
    })
  )
  cor_sig <- cor_dt[toupper(gene) %in% signature_genes]
  
  # ---------------------------
  # Merge DE + correlation
  # ---------------------------
  merged <- merge(de_sig, cor_sig, by = "gene", all = FALSE)
  merged[, direction_match := sign(logFC) == sign(rho)]
  
  # ---------------------------
  # Summary
  # ---------------------------
  summary_dt <- data.table(
    signature = signature_name,
    cohort = cohort_name,
    n_signature_genes_input = length(signature_genes),
    n_signature_genes_detected = nrow(merged),
    n_DE_FDR_lt_0.05 = sum(merged$adj.P.Val < 0.05, na.rm = TRUE),
    n_DE_FDR_lt_0.05_and_positive = sum(merged$adj.P.Val < 0.05 & merged$logFC > 0, na.rm = TRUE),
    n_direction_match = sum(merged$direction_match, na.rm = TRUE),
    median_rho = median(merged$rho, na.rm = TRUE)
  )
  
  # ---------------------------
  # Plot
  # ---------------------------
  p <- ggplot(merged, aes(x = logFC, y = rho)) +
    geom_point(shape = 1, size = 2.2, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(cohort_name, ": ", signature_name),
      x = "log2FC (ERBB2-high vs ERBB2-low)",
      y = "Spearman rho with ERBB2 abundance"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # ---------------------------
  # Save outputs
  # ---------------------------
  prefix <- paste0(cohort_name, "_", signature_name)
  
  fwrite(de_full, file.path(out_dir, paste0(prefix, "_DE_full.csv")))
  fwrite(de_sig, file.path(out_dir, paste0(prefix, "_DE_signature_only.csv")))
  fwrite(cor_sig, file.path(out_dir, paste0(prefix, "_ERBB2_correlation_signature_only.csv")))
  fwrite(merged, file.path(out_dir, paste0(prefix, "_DE_COR_combined.csv")))
  fwrite(summary_dt, file.path(out_dir, paste0(prefix, "_summary.csv")))
  
  ggsave(
    filename = file.path(out_dir, paste0(prefix, "_DE_vs_COR.png")),
    plot = p, width = 6.5, height = 5.2, dpi = 300
  )
  
  print(summary_dt)
  invisible(summary_dt)
}

# -------------------------------
# LOAD MATRICES ONCE
# -------------------------------
ov_mat   <- prepare_matrix(ov_fp, "OVARIAN")
ucec_mat <- prepare_matrix(ucec_fp, "UCEC")

# -------------------------------
# LOAD SIGNATURES
# -------------------------------
sig_list <- lapply(sig_files, load_signature)

for (nm in names(sig_list)) {
  cat("Loaded signature", nm, ":", length(sig_list[[nm]]), "genes\n")
}

# -------------------------------
# RUN ALL COMBINATIONS
# -------------------------------
all_summaries <- list()

for (sig_name in names(sig_list)) {
  sig_genes <- sig_list[[sig_name]]
  
  all_summaries[[paste0("OVARIAN_", sig_name)]] <-
    run_projection_analysis(ov_mat, "OVARIAN", sig_genes, sig_name, out_dir)
  
  all_summaries[[paste0("UCEC_", sig_name)]] <-
    run_projection_analysis(ucec_mat, "UCEC", sig_genes, sig_name, out_dir)
}

final_summary <- rbindlist(all_summaries, fill = TRUE)
fwrite(final_summary, file.path(out_dir, "FINAL_SUMMARY_all_signatures.csv"))

cat("\nDONE.\nOutputs written to:\n", out_dir, "\n")
print(final_summary)