suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

# ============================================================
# Summarize ProteoMixture results across proteomic cohorts
# ============================================================

# -----------------------------
# User paths: update these if needed
# -----------------------------
files <- list(
  RajKumar = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/raj kumar et al/aim2_outputs/RajKumar_PROT_for_ProteoMixture_ProteoMixture_ssGSEA_results.csv",
  Krug     = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/krug et al/aim2_outputs/Krug_PROT_for_ProteoMixture_ProteoMixture_ssGSEA_results.csv",
  Mertins  = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/mertins et al/aim2_outputs/Mertins_PROT_for_ProteoMixture_ProteoMixture_ssGSEA_results.csv"
)

out_dir <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/proteomixture_summary_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper function
# -----------------------------
analyze_proteomixture <- function(fp, cohort_name, out_dir) {
  dt <- fread(fp)
  
  # Standardize expected columns
  needed <- c("Tumor", "Immune", "Stroma")
  if (!all(needed %in% names(dt))) {
    stop(
      cohort_name, ": missing expected columns. Found: ",
      paste(names(dt), collapse = ", ")
    )
  }
  
  # sample names column, if present
  sample_col <- NULL
  if ("sample_names" %in% names(dt)) sample_col <- "sample_names"
  if (is.null(sample_col) && "V1" %in% names(dt)) sample_col <- "V1"
  if (is.null(sample_col)) {
    dt[, sample_names := paste0(cohort_name, "_", seq_len(.N))]
    sample_col <- "sample_names"
  }
  
  # Correlations
  rho_tumor_immune  <- suppressWarnings(cor(dt$Tumor, dt$Immune, method = "spearman", use = "complete.obs"))
  rho_tumor_stroma  <- suppressWarnings(cor(dt$Tumor, dt$Stroma, method = "spearman", use = "complete.obs"))
  
  cor_summary <- data.table(
    cohort = cohort_name,
    n_samples = nrow(dt),
    rho_tumor_vs_immune = rho_tumor_immune,
    rho_tumor_vs_stroma = rho_tumor_stroma
  )
  
  # Scatter: Tumor vs Immune
  p1 <- ggplot(dt, aes(x = Immune, y = Tumor)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(cohort_name, ": Tumor vs Immune"),
      subtitle = paste0("Spearman \u03c1 = ", round(rho_tumor_immune, 3)),
      x = "Immune score",
      y = "Tumor score"
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # Scatter: Tumor vs Stroma
  p2 <- ggplot(dt, aes(x = Stroma, y = Tumor)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(cohort_name, ": Tumor vs Stroma"),
      subtitle = paste0("Spearman \u03c1 = ", round(rho_tumor_stroma, 3)),
      x = "Stroma score",
      y = "Tumor score"
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # Save individual plots
  ggsave(
    file.path(out_dir, paste0(cohort_name, "_Tumor_vs_Immune.png")),
    p1, width = 6.5, height = 5.2, dpi = 300
  )
  ggsave(
    file.path(out_dir, paste0(cohort_name, "_Tumor_vs_Stroma.png")),
    p2, width = 6.5, height = 5.2, dpi = 300
  )
  
  # Save one combined panel per cohort
  png(file.path(out_dir, paste0(cohort_name, "_ProteoMixture_scatter_panel.png")),
      width = 1800, height = 800, res = 150)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
  
  list(
    summary = cor_summary,
    p_tumor_immune = p1,
    p_tumor_stroma = p2
  )
}

# -----------------------------
# Run all cohorts
# -----------------------------
results <- list()
summary_list <- list()

for (nm in names(files)) {
  message("Running: ", nm)
  results[[nm]] <- analyze_proteomixture(files[[nm]], nm, out_dir)
  summary_list[[nm]] <- results[[nm]]$summary
}

cor_table <- rbindlist(summary_list)
fwrite(cor_table, file.path(out_dir, "ProteoMixture_correlation_summary_all_cohorts.tsv"), sep = "\t")

# -----------------------------
# Combined summary figure
# -----------------------------
cor_long <- melt(
  cor_table,
  id.vars = c("cohort", "n_samples"),
  measure.vars = c("rho_tumor_vs_immune", "rho_tumor_vs_stroma"),
  variable.name = "comparison",
  value.name = "spearman_rho"
)

cor_long[, comparison := fifelse(
  comparison == "rho_tumor_vs_immune",
  "Tumor vs Immune",
  "Tumor vs Stroma"
)]

p_summary <- ggplot(cor_long, aes(x = cohort, y = spearman_rho)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  facet_wrap(~comparison, nrow = 1) +
  theme_classic(base_size = 14) +
  labs(
    title = "ProteoMixture correlation summary across proteomic cohorts",
    x = "Cohort",
    y = "Spearman correlation"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  file.path(out_dir, "ProteoMixture_correlation_summary_all_cohorts.png"),
  p_summary, width = 10, height = 4.8, dpi = 300
)

message("Done. Outputs written to: ", out_dir)
print(cor_table)