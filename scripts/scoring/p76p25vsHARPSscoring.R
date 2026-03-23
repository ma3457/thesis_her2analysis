suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------------- paths ----------------
base_out <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs"

p76p25_fp <- file.path(base_out, "Wolf_P76P25_scoring", "Wolf_P76P25_scores_meanZ.tsv")
harps_fp  <- file.path(base_out, "Wolf_HArPS_score_RNA", "Wolf_HARPS_score_per_sample.csv")

stopifnot(file.exists(p76p25_fp), file.exists(harps_fp))

# ---------------- read ----------------
p76p25 <- fread(p76p25_fp)
harps  <- fread(harps_fp)

# Make join key consistent (this fixes your double vs character join error)
# Adjust the column names below if yours differ.
# P76/P25 file: usually has sample_id + P76_meanZ + P25_meanZ (or similar)
# HARPS file:   usually has sample_id + HARPS_score (or similar)

# If sample_id column name differs, rename it here:
if (!"sample_id" %in% names(p76p25)) {
  # common alternatives
  alt <- intersect(names(p76p25), c("Patient Identifier","patient_id","id","sample","Sample","SampleID"))
  if (length(alt) == 0) stop("Can't find sample_id in p76p25. Columns: ", paste(names(p76p25), collapse=", "))
  setnames(p76p25, alt[1], "sample_id")
}
if (!"sample_id" %in% names(harps)) {
  alt <- intersect(names(harps), c("Patient Identifier","patient_id","id","sample","Sample","SampleID"))
  if (length(alt) == 0) stop("Can't find sample_id in harps. Columns: ", paste(names(harps), collapse=", "))
  setnames(harps, alt[1], "sample_id")
}

p76p25[, sample_id := as.character(sample_id)]
harps[,  sample_id := as.character(sample_id)]

# Identify score columns robustly
pick_col <- function(dt, patterns) {
  hits <- unlist(lapply(patterns, function(p) grep(p, names(dt), value = TRUE, ignore.case = TRUE)))
  hits <- unique(hits)
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

p76_col  <- pick_col(p76p25, c("^P76", "P76_mean", "P76.*score", "P76.*meanZ"))
p25_col  <- pick_col(p76p25, c("^P25", "P25_mean", "P25.*score", "P25.*meanZ"))
harp_col <- pick_col(harps,  c("HARPS_score", "HARP.*score", "score"))

stopifnot(!is.na(p76_col), !is.na(p25_col), !is.na(harp_col))

setnames(p76p25, c(p76_col, p25_col), c("P76_score", "P25_score"))
setnames(harps,  harp_col, "HARPS_score")

# ---------------- merge ----------------
dt <- merge(p76p25[, .(sample_id, P76_score, P25_score)],
            harps[,  .(sample_id, HARPS_score)],
            by = "sample_id")

cat("Merged N =", nrow(dt), "\n")

# ---------------- stats ----------------
rho_p76 <- suppressWarnings(cor(dt$P76_score, dt$HARPS_score, method="spearman", use="complete.obs"))
rho_p25 <- suppressWarnings(cor(dt$P25_score, dt$HARPS_score, method="spearman", use="complete.obs"))

cat("Spearman rho:\n")
cat("  P76 vs HARPS =", round(rho_p76, 3), "\n")
cat("  P25 vs HARPS =", round(rho_p25, 3), "\n")

# ---------------- plots ----------------
out_dir <- file.path(base_out, "Wolf_P76P25_vs_HARPS")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

make_scatter <- function(x, y, xlab, ylab, title, out_png) {
  p <- ggplot(dt, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(shape = 1, size = 2.2, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    theme_classic(base_size = 16) +
    labs(x = xlab, y = ylab, title = title) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(out_png, p, width = 7.6, height = 5.2, dpi = 300)
}

make_scatter("P76_score", "HARPS_score",
             "P76 pathway score (mean z)", "HARPS-score",
             sprintf("Wolf: P76 vs HARPS-score (rho=%.3f)", rho_p76),
             file.path(out_dir, "Wolf_P76_vs_HARPS_scatter.png"))

make_scatter("P25_score", "HARPS_score",
             "P25 pathway score (mean z)", "HARPS-score",
             sprintf("Wolf: P25 vs HARPS-score (rho=%.3f)", rho_p25),
             file.path(out_dir, "Wolf_P25_vs_HARPS_scatter.png"))

# Optional: one combined panel-style plot (P76 on x, P25 on y colored by HARPS)
p_combo <- ggplot(dt, aes(x = P76_score, y = P25_score, color = HARPS_score)) +
  geom_point(size = 2.2, alpha = 0.9) +
  theme_classic(base_size = 16) +
  labs(title = "Wolf: P76 vs P25 (colored by HARPS-score)",
       x = "P76 pathway score (mean z)",
       y = "P25 pathway score (mean z)")
ggsave(file.path(out_dir, "Wolf_P76_vs_P25_colored_by_HARPS.png"),
       p_combo, width = 7.6, height = 5.2, dpi = 300)

# Save merged table for slides
fwrite(dt, file.path(out_dir, "Wolf_P76P25_vs_HARPS_merged.tsv"), sep = "\t")

cat("✅ Wrote outputs to:\n", out_dir, "\n")
