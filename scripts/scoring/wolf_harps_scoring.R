# ============================ wolf_harps_scoring.R ============================
# Wolf (I-SPY2 / GSE194040): make HArPS labels from RPS-5, run RNA DE
# (HARP+ vs HARP- within HER2-), then build a sample-level "HARPS-score".
#
# Outputs -> <script_dir>/aim2_outputs/Wolf_HARPS_score_RNA/
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(limma)
  library(stringr)
  library(ggplot2)
})

# ---------------- script directory ----------------
script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) && file.exists(f)) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(".", winslash = "/")
})

# ==============================================================================
# CONFIG
# ==============================================================================

SUPP3_XLSX <- file.path(script_dir, "NIHMS1829047-supplement-3.xlsx")
RNA_FILE   <- file.path(script_dir, "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt")

OUT_ROOT <- file.path(script_dir, "aim2_outputs")
OUT_DIR  <- file.path(OUT_ROOT, "Wolf_HARPS_score_RNA")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Use this threshold to define the scoring gene set
FDR_MAX <- 0.05

# If LFC1 is too strict, keep this at 0 (or try 0.25 / 0.5 later)
LFC_MIN_FOR_SET <- 0

# ==============================================================================
# Helpers
# ==============================================================================

norm_rps5 <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[\u2212\u2013\u2014]", "-", x)  # unicode dashes -> minus
  x <- gsub("\\s+", "", x)
  x <- gsub("_", "/", x)
  x
}

her2_status_from_rps5 <- function(rps5_raw) {
  s <- norm_rps5(rps5_raw)
  out <- rep(NA_character_, length(s))
  out[grepl("^her2\\+/", s)] <- "HER2+"
  out[grepl("^her2-/",  s)]  <- "HER2-"
  factor(out, levels = c("HER2-","HER2+"))
}

map_harp_from_rps5 <- function(rps5_raw) {
  s <- norm_rps5(rps5_raw)
  out <- rep(NA_character_, length(s))
  
  is_her2_pos <- grepl("^her2\\+/", s)
  is_her2_neg <- grepl("^her2-/",  s)
  imm_pos     <- grepl("immune\\+", s)
  imm_neg     <- grepl("immune\\-", s)
  drd_pos     <- grepl("drd(\\.v\\d+)?\\+", s)
  drd_neg     <- grepl("drd(\\.v\\d+)?\\-", s)
  no_drd      <- !grepl("drd", s)
  
  # HArPS+ : HER2− and (immune+ OR DRD+)
  out[is_her2_neg & (imm_pos | drd_pos)] <- "HARP+"
  # HArPS− : HER2− and immune− and (DRD− OR no DRD annotation)
  out[is_her2_neg & imm_neg & (drd_neg | no_drd)] <- "HARP-"
  
  # Exclude HER2+ and missing RPS-5 from HArPS assignment
  out[is.na(rps5_raw) | is_her2_pos] <- NA_character_
  
  factor(out, levels = c("HARP-","HARP+"))
}

# Wide reader that tolerates messy headers (first col = gene IDs, rest = samples)
read_expr_wide <- function(path, feature_name = "gene_symbol") {
  stopifnot(file.exists(path))
  
  hdr_line   <- tryCatch(readLines(path, n = 1), error = function(e) character(0))
  hdr_tokens <- if (length(hdr_line)) strsplit(hdr_line, "\t", fixed = TRUE)[[1]] else character(0)
  
  dt <- fread(path, sep = "\t", header = FALSE, fill = TRUE, quote = "", check.names = FALSE)
  
  if (length(hdr_tokens) == ncol(dt)) {
    setnames(dt, hdr_tokens)
  } else if (length(hdr_tokens) + 1L == ncol(dt)) {
    setnames(dt, c(feature_name, hdr_tokens))
  } else {
    setnames(dt, c(feature_name, paste0("Sample_", seq_len(ncol(dt) - 1L))))
  }
  
  id_col <- names(dt)[1]
  if (is.na(id_col) || id_col == "") { setnames(dt, 1, feature_name); id_col <- feature_name }
  if (!id_col %in% c(feature_name,"gene_symbol","Gene","feature_id","ID","Feature")) {
    setnames(dt, 1, feature_name); id_col <- feature_name
  }
  
  dt <- dt[!is.na(get(id_col)) & trimws(get(id_col)) != "", ]
  dt <- dt[!duplicated(get(id_col)), ]
  
  samp_cols <- setdiff(names(dt), id_col)
  if (length(samp_cols) == 0) stop("No sample columns found after header repair.")
  
  # safer numeric coercion
  for (nm in samp_cols) {
    if (!is.numeric(dt[[nm]])) {
      v <- as.character(dt[[nm]])
      v <- trimws(v)
      v[v %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
      suppressWarnings(dt[[nm]] <- as.numeric(v))
    }
  }
  
  keep_rows <- dt[, rowSums(as.data.frame(lapply(.SD, function(z) !is.na(z)))) > 0, .SDcols = samp_cols]
  if (any(!keep_rows)) dt <- dt[keep_rows]
  
  E  <- as.matrix(dt[, ..samp_cols])
  rn <- make.unique(as.character(dt[[id_col]]))
  n  <- min(nrow(E), length(rn))
  E  <- E[seq_len(n), , drop = FALSE]
  rownames(E) <- rn[seq_len(n)]
  storage.mode(E) <- "numeric"
  
  # IMPORTANT: keep numeric-looking patient IDs as strings
  colnames(E) <- as.character(colnames(E))
  
  E
}

run_limma_two_group <- function(expr, meta_dt, group_col, pos_label = "HARP+", neg_label = "HARP-") {
  # ensure character IDs everywhere (prevents bmerge join-type issues later too)
  meta_dt[, sample_id := as.character(sample_id)]
  colnames(expr) <- as.character(colnames(expr))
  
  keep <- intersect(colnames(expr), meta_dt$sample_id)
  if (!length(keep)) stop("No overlap between expression columns and metadata sample_id.")
  
  X   <- expr[, keep, drop = FALSE]
  met <- meta_dt[match(colnames(X), meta_dt$sample_id)]
  grp <- factor(met[[group_col]], levels = c(neg_label, pos_label))
  
  vals <- as.numeric(X[is.finite(X)])
  if (length(vals) && (quantile(vals, 0.99, na.rm = TRUE) > 30 || max(vals, na.rm = TRUE) > 1000)) {
    X <- log2(X + 1)
  }
  
  zv <- apply(X, 1, sd, na.rm = TRUE) == 0
  if (any(zv)) X <- X[!zv, , drop = FALSE]
  
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- c("neg","pos")
  
  fit  <- lmFit(X, design)
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(pos - neg, levels = design)),
                 trend = TRUE, robust = TRUE)
  
  tt <- topTable(fit2, number = Inf, sort.by = "P")
  data.table(
    feature_id = rownames(tt),
    logFC      = tt$logFC,
    AveExpr    = tt$AveExpr,
    t          = tt$t,
    P.Value    = tt$P.Value,
    adj.P.Val  = tt$adj.P.Val
  )
}

# Score = mean(z(up)) - mean(z(down)), where up/down are based on DE direction
harps_score_from_de <- function(expr_mat, de_dt, fdr_max = 0.05, lfc_min = 0) {
  stopifnot(is.matrix(expr_mat))
  stopifnot(all(c("feature_id","logFC","adj.P.Val") %in% names(de_dt)))
  
  hits <- de_dt[!is.na(adj.P.Val) & adj.P.Val < fdr_max & abs(logFC) >= lfc_min]
  
  up   <- hits[logFC > 0, unique(feature_id)]
  down <- hits[logFC < 0, unique(feature_id)]
  
  up   <- intersect(up, rownames(expr_mat))
  down <- intersect(down, rownames(expr_mat))
  
  # gene-wise z across samples
  Z <- t(scale(t(expr_mat)))
  
  z_up   <- if (length(up)   >= 1) colMeans(Z[up, , drop = FALSE],   na.rm = TRUE) else rep(0, ncol(Z))
  z_down <- if (length(down) >= 1) colMeans(Z[down, , drop = FALSE], na.rm = TRUE) else rep(0, ncol(Z))
  
  list(score = z_up - z_down, up = up, down = down)
}

save_box <- function(dt, x, y, title, out_png) {
  p <- ggplot(dt, aes(x = .data[[x]], y = .data[[y]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.18, height = 0, alpha = 0.7, size = 1.4) +
    theme_classic(base_size = 14) +
    labs(title = title, x = NULL, y = y) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(out_png, p, width = 7.2, height = 5.2, dpi = 300)
}

# ==============================================================================
# STEP 1: Labels from Supplement 3 (RPS-5)
# ==============================================================================

stopifnot(file.exists(SUPP3_XLSX))
supp3 <- as.data.table(read_xlsx(SUPP3_XLSX, sheet = 1))

meta <- supp3[, .(`Patient Identifier`, `RPS-5`)]
setnames(meta, c("Patient Identifier","RPS-5"), c("sample_id","RPS5"))

# force IDs to character so merges never break
meta[, sample_id := as.character(sample_id)]

meta[, HER2_status := her2_status_from_rps5(RPS5)]
meta[, HARP        := map_harp_from_rps5(RPS5)]

meta_labels <- meta[, .(sample_id, HER2_status, HARP, RPS5)]
fwrite(meta_labels, file.path(OUT_DIR, "labels_from_RPS5.csv"))

meta_harp <- meta[!is.na(HARP), .(sample_id, HARP, RPS5)]
meta_harp[, sample_id := as.character(sample_id)]
fwrite(meta_harp, file.path(OUT_DIR, "metadata_HARP_only.csv"))

cat("HArPS group counts (HER2- only):\n")
print(table(meta_harp$HARP, useNA = "ifany"))

# ==============================================================================
# STEP 2: RNA DE (HARP+ vs HARP- within HER2-)
# ==============================================================================

stopifnot(file.exists(RNA_FILE))
E_rna <- read_expr_wide(RNA_FILE, feature_name = "gene_symbol")

cat("\n--- Running limma: HARP+ vs HARP- (HER2- only) ---\n")
rna_tt <- run_limma_two_group(E_rna, meta_harp, group_col = "HARP", pos_label = "HARP+", neg_label = "HARP-")

fwrite(rna_tt, file.path(OUT_DIR, "LIMMA_all.csv"))
fwrite(rna_tt[adj.P.Val < FDR_MAX], file.path(OUT_DIR, "LIMMA_FDR0.05.csv"))

cat("RNA DE rows:\n")
cat("  all = ", nrow(rna_tt), "\n", sep = "")
cat("  FDR<0.05 = ", nrow(rna_tt[adj.P.Val < FDR_MAX]), "\n", sep = "")

# ==============================================================================
# STEP 3: Build HARPS-score from the DE result
# ==============================================================================

sc <- harps_score_from_de(
  expr_mat = E_rna,
  de_dt    = rna_tt,
  fdr_max  = FDR_MAX,
  lfc_min  = LFC_MIN_FOR_SET
)

score_dt <- data.table(
  sample_id   = as.character(colnames(E_rna)),
  HARPS_score = as.numeric(sc$score)
)

meta_labels[, sample_id := as.character(sample_id)]
score_dt <- merge(score_dt, meta_labels[, .(sample_id, HER2_status, HARP)], by = "sample_id", all.x = TRUE)

fwrite(score_dt, file.path(OUT_DIR, "Wolf_HARPS_score_per_sample.csv"))
fwrite(data.table(gene = sc$up),   file.path(OUT_DIR, "HARPS_score_genes_UP.csv"))
fwrite(data.table(gene = sc$down), file.path(OUT_DIR, "HARPS_score_genes_DOWN.csv"))

cat("\nHARPS-score gene set sizes:\n")
cat("  UP = ", length(sc$up), "\n", sep = "")
cat("  DOWN = ", length(sc$down), "\n", sep = "")

# ==============================================================================
# STEP 4: Plots
# ==============================================================================

plot_dt1 <- score_dt[!is.na(HARP)]
if (nrow(plot_dt1)) {
  save_box(
    dt = plot_dt1,
    x = "HARP",
    y = "HARPS_score",
    title = "Wolf (HER2- only): HARPS-score by HArPS label",
    out_png = file.path(OUT_DIR, "HARPS_score_by_HARP.png")
  )
}

plot_dt2 <- score_dt[!is.na(HER2_status)]
if (nrow(plot_dt2)) {
  save_box(
    dt = plot_dt2,
    x = "HER2_status",
    y = "HARPS_score",
    title = "Wolf: HARPS-score by HER2 status (from RPS-5)",
    out_png = file.path(OUT_DIR, "HARPS_score_by_HER2status.png")
  )
}

if ("ERBB2" %in% rownames(E_rna)) {
  erbb2_dt <- data.table(
    sample_id = as.character(colnames(E_rna)),
    ERBB2     = as.numeric(E_rna["ERBB2", ])
  )
  erbb2_dt <- merge(erbb2_dt, score_dt, by = "sample_id", all.x = TRUE)
  
  p <- ggplot(erbb2_dt, aes(x = ERBB2, y = HARPS_score)) +
    geom_point(shape = 1, size = 2.0, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    theme_classic(base_size = 14) +
    labs(title = "Wolf: ERBB2 expression vs HARPS-score", x = "ERBB2 expression", y = "HARPS-score") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(OUT_DIR, "HARPS_score_vs_ERBB2.png"), p, width = 7.2, height = 5.2, dpi = 300)
}

# ==============================================================================
# RUN LOG
# ==============================================================================

log_lines <- c(
  "RUN LOG — Wolf HARPS-score",
  paste0("script_dir: ", script_dir),
  paste0("SUPP3_XLSX: ", basename(SUPP3_XLSX)),
  paste0("RNA_FILE: ", basename(RNA_FILE)),
  paste0("OUT_DIR: ", OUT_DIR),
  paste0("FDR_MAX: ", FDR_MAX),
  paste0("LFC_MIN_FOR_SET: ", LFC_MIN_FOR_SET),
  paste0("DE rows (all): ", nrow(rna_tt)),
  paste0("DE rows (FDR<0.05): ", nrow(rna_tt[adj.P.Val < FDR_MAX])),
  paste0("HARPS UP genes used: ", length(sc$up)),
  paste0("HARPS DOWN genes used: ", length(sc$down))
)

writeLines(log_lines, file.path(OUT_DIR, "RUN_LOG.txt"))
cat("\n✅ Done. Outputs in:\n", OUT_DIR, "\n", sep = "")
