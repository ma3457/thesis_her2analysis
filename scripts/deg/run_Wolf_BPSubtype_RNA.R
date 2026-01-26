# run_Wolf_BPsubtype_RNA.R
# Wolf (GSE194040 / I-SPY2) RNA expression: HER2pos vs HER2neg using BP-subtype labels.
#
# Outputs (PNG-only):
#   - aim2_outputs/<tag>/volcano.png
#   - aim2_outputs/<tag>/MA.png
#   - aim2_outputs/<tag>/qc_ERBB2_boxplot.png
#   - aim2_outputs/<tag>/qc_GRB7_boxplot.png
#   - aim2_outputs/<tag>/CORR_HER2pos_vs_HER2neg_filtered.png (if enough genes)
#   - aim2_outputs/<tag>/RUN_LOG.txt
# Tables:
#   - ALL_HER2pos_vs_HER2neg_DE.tsv
#   - FILT_HER2pos_vs_HER2neg_DE_FDR0.05.tsv
#   - FILT_HER2pos_vs_HER2neg_DE_FDR0.05_LFC1.tsv
# Meta artifacts:
#   - <tag>_RNA_meta_template.csv (template)
#   - aim2_outputs/<tag>/<tag>_RNA_meta_USED.tsv
#   - aim2_outputs/<tag>/<tag>_RNA_contrast_USED.txt

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(pheatmap)
  library(stringr)
})

# --- discover the folder this script lives in (no setwd required) ---
script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) && file.exists(f)) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(".", winslash = "/")
})

# --------- paths ----------
expr_file <- file.path(
  script_dir,
  "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt"
)
lab_file  <- file.path(script_dir, "labels_from_BP_subtype.csv")

out_tag <- "GSE194040_ISPY2_mRNA_BPsubtype"
out_dir <- file.path(script_dir, "aim2_outputs", out_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------- meta artifacts ----------
meta_template_csv <- file.path(script_dir, paste0(out_tag, "_RNA_meta_template.csv"))
meta_used_tsv     <- file.path(out_dir,  paste0(out_tag, "_RNA_meta_USED.tsv"))
contrast_used_txt <- file.path(out_dir,  paste0(out_tag, "_RNA_contrast_USED.txt"))

wrap_title <- function(s, width = 60) paste(strwrap(s, width = width), collapse = "\n")
plot_with_margins <- function(expr) {
  op <- par(mar = c(5, 5, 5, 2))
  on.exit(par(op), add = TRUE)
  force(expr)
}

stopifnot(file.exists(expr_file), file.exists(lab_file))

# If a prior run crashed mid-plot, a device may still be open
graphics.off()

# ---------------- read expression ----------------
x <- fread(expr_file, check.names = FALSE)

# enforce a clean gene ID column
if (is.na(names(x)[1]) || names(x)[1] == "") setnames(x, 1, "gene_symbol")
if (!names(x)[1] %in% c("gene_symbol", "Gene", "feature_id", "ID", "Feature")) {
  setnames(x, 1, "gene_symbol")
}
x <- x[!duplicated(x$gene_symbol), ]

# coerce numeric sample columns safely
samp_cols <- setdiff(names(x), "gene_symbol")
for (nm in samp_cols) {
  if (!is.numeric(x[[nm]])) {
    v <- as.character(x[[nm]])
    v[is.na(v)] <- ""
    suppressWarnings(num <- as.numeric(v))
    if (sum(nchar(v) > 0 & is.na(num)) == 0) x[[nm]] <- num
  }
}

E <- as.matrix(x[, ..samp_cols])
rownames(E) <- x$gene_symbol
storage.mode(E) <- "numeric"

# ---------------- read BP-subtype labels ----------------
lab <- fread(lab_file, check.names = FALSE)
stopifnot(all(c("sample_id", "group") %in% names(lab)))

lab$sample_id <- str_trim(as.character(lab$sample_id))
lab$group     <- str_trim(as.character(lab$group))

# write a simple meta template (all labeled samples). Optional to edit.
if (!file.exists(meta_template_csv)) {
  meta_template <- data.table(
    SampleID  = lab$sample_id,
    Group_raw = lab$group,  # Positive / Negative as provided
    Notes     = ""
  )
  fwrite(meta_template, meta_template_csv)
  message("Wrote meta template: ", meta_template_csv)
} else {
  message("Meta template exists: ", meta_template_csv)
}

# ---------------- align by exact ID; fallback to shared numeric key ----------------
keep <- intersect(colnames(E), lab$sample_id)

if (!length(keep)) {
  key_lab  <- str_extract(lab$sample_id, "\\d{5,}")
  key_expr <- str_extract(colnames(E),   "\\d{5,}")
  
  m <- merge(
    data.table(sample_id = lab$sample_id, key = key_lab)[!is.na(key)],
    data.table(expr_col  = colnames(E),   key = key_expr)[!is.na(key)],
    by = "key"
  )
  
  if (nrow(m)) {
    m <- m[, .SD[.N == 1], by = key]
    lab <- merge(lab, unique(m[, .(sample_id, expr_col)]), by = "sample_id", all.x = TRUE)
    lab$sample_id <- ifelse(!is.na(lab$expr_col), lab$expr_col, lab$sample_id)
    lab$expr_col  <- NULL
  }
  
  fwrite(m, file.path(out_dir, "mapping_help_labels_to_expr.tsv"), sep = "\t")
  keep <- intersect(colnames(E), lab$sample_id)
}

if (!length(keep)) stop("No overlap between expression columns and labels$sample_id")

# finalize alignment + drop NA/other labels
E   <- E[, keep, drop = FALSE]
lab <- lab[match(colnames(E), lab$sample_id), , drop = FALSE]

ok  <- !is.na(lab$group) & lab$group %in% c("Positive", "Negative")
E   <- E[, ok, drop = FALSE]
lab <- lab[ok, , drop = FALSE]

# groups
grp <- factor(ifelse(lab$group == "Positive", "HER2pos", "HER2neg"),
              levels = c("HER2neg", "HER2pos"))

n_pos <- sum(grp == "HER2pos")
n_neg <- sum(grp == "HER2neg")
if (n_pos < 3 || n_neg < 3) stop("Groups too small (need ≥3 per arm).")

# write meta USED + contrast USED
meta_used <- data.table(
  SampleID  = colnames(E),
  Group     = as.character(grp),
  Group_raw = lab$group
)
fwrite(meta_used, meta_used_tsv, sep = "\t")

writeLines(c(
  "Mode: HER2",
  "Contrast: HER2pos - HER2neg",
  paste0("Cohort: ", out_tag),
  paste0("n_pos: ", n_pos),
  paste0("n_neg: ", n_neg)
), contrast_used_txt)

message("Wrote meta USED: ", meta_used_tsv)
message("Wrote contrast USED: ", contrast_used_txt)

# quick log2 if values look raw-ish
vals <- as.numeric(E[is.finite(E)])
if (quantile(vals, 0.99, na.rm = TRUE) > 30 || max(vals, na.rm = TRUE) > 1000) {
  E <- log2(E + 1)
}

# drop zero-variance
zv <- apply(E, 1, sd, na.rm = TRUE) == 0
if (any(zv)) E <- E[!zv, , drop = FALSE]

# ---------------- LIMMA (HER2pos - HER2neg) ----------------
design <- model.matrix(~ 0 + grp)
colnames(design) <- c("HER2neg", "HER2pos")

fit  <- lmFit(E, design)
fit2 <- eBayes(
  contrasts.fit(fit, makeContrasts(HER2pos - HER2neg, levels = design)),
  trend  = TRUE,
  robust = TRUE
)

tt <- topTable(fit2, number = Inf, sort.by = "P")
tt$feature_id <- rownames(tt)

# significance flags and counts
sig_flag <- tt$adj.P.Val < 0.05
n_fdr    <- sum(sig_flag, na.rm = TRUE)

fc2_flag <- sig_flag & abs(tt$logFC) >= 1
tab_fc2  <- tt[fc2_flag, ]
n_fc2    <- nrow(tab_fc2)

# ---------------- outputs ----------------
fwrite(tt, file.path(out_dir, "ALL_HER2pos_vs_HER2neg_DE.tsv"), sep = "\t")

tab_fdr <- tt[sig_flag, c("feature_id", "logFC", "adj.P.Val")]
fwrite(tab_fdr, file.path(out_dir, "FILT_HER2pos_vs_HER2neg_DE_FDR0.05.tsv"), sep = "\t")

fwrite(tab_fc2[, c("feature_id", "logFC", "adj.P.Val")],
       file.path(out_dir, "FILT_HER2pos_vs_HER2neg_DE_FDR0.05_LFC1.tsv"), sep = "\t")

# correlation heatmap from FC2-filtered genes (needs ≥3 genes)
if (nrow(tab_fc2) >= 3) {
  keep_genes <- intersect(tab_fc2$feature_id, rownames(E))
  if (length(keep_genes) >= 3) {
    M  <- E[keep_genes, , drop = FALSE]
    cm <- cor(t(M), use = "pairwise.complete.obs")
    png(file.path(out_dir, "CORR_HER2pos_vs_HER2neg_filtered.png"),
        width = 1400, height = 1200, res = 180)
    pheatmap(cm, main = "Correlation — HER2pos vs HER2neg (FDR<0.05 & |log2FC|>=1)")
    dev.off()
  }
}

# ---------------- QC boxplots (PNG-only) ----------------
for (g in c("ERBB2", "GRB7")) {
  if (g %in% rownames(E)) {
    ttl <- wrap_title(sprintf("%s %s | BP-subtype (n+=%d, n-=%d)", out_tag, g, n_pos, n_neg), 55)
    fn  <- file.path(out_dir, sprintf("qc_%s_boxplot.png", g))
    
    png(fn, width = 1400, height = 1100, res = 150)
    plot_with_margins({
      boxplot(split(as.numeric(E[g, ]), grp),
              main = ttl, ylab = "Expression", xlab = "HER2 group",
              col = "grey90", border = "grey30", cex.main = 0.95)
    })
    dev.off()
  }
}

# ---------------- Volcano (PNG-only; -log10 FDR) ----------------
plot_volcano <- function() {
  ttl <- wrap_title(sprintf("%s Volcano (BP-subtype, n+=%d, n-=%d)", out_tag, n_pos, n_neg), 55)
  y <- -log10(pmax(tt$adj.P.Val, 1e-300))
  
  plot(tt$logFC, y,
       xlab = "log2FC (pos vs neg)",
       ylab = expression(-log[10]("FDR")),
       main = ttl,
       pch = 21,
       bg  = ifelse(sig_flag, "black", "white"),
       col = "black",
       cex.main = 0.95)
  
  abline(v = c(-1, 1), lty = 2, col = "grey75")
  abline(h = -log10(0.05), lty = 2, col = "grey75")
  
  tops <- unique(c(head(tt[order(tt$adj.P.Val), "feature_id"], 5), "ERBB2", "GRB7"))
  tops <- tops[tops %in% tt$feature_id]
  sel  <- tt$feature_id %in% tops
  if (any(sel)) {
    text(tt$logFC[sel], y[sel], labels = tt$feature_id[sel], pos = 4, cex = 0.8)
  }
}

png(file.path(out_dir, "volcano.png"), width = 1400, height = 1100, res = 150)
plot_with_margins(plot_volcano())
dev.off()

# ---------------- MA (PNG-only) ----------------
plot_ma <- function() {
  ttl <- wrap_title(sprintf("%s MA (BP-subtype, n+=%d, n-=%d)", out_tag, n_pos, n_neg), 55)
  plot(fit$Amean, fit2$coefficients[, 1],
       xlab = "AveExpr", ylab = "log2FC",
       main = ttl,
       pch = 21, bg = "white", col = "black",
       cex.main = 0.95)
  abline(h = 0, lty = 2, col = "grey50")
}

png(file.path(out_dir, "MA.png"), width = 1400, height = 1100, res = 150)
plot_with_margins(plot_ma())
dev.off()

# ---------------- run log ----------------
pick_line <- function(df, g) {
  i <- which(df$feature_id == g)[1]
  if (!is.na(i)) sprintf("%s: log2FC=%.3f  FDR=%.3g", g, df$logFC[i], df$adj.P.Val[i]) else sprintf("%s: NA", g)
}

log_txt <- c(
  sprintf("RUN LOG — %s", out_tag),
  sprintf("Matrix: %s", basename(expr_file)),
  sprintf("Labels: %s", basename(lab_file)),
  sprintf("Samples kept: %d", ncol(E)),
  sprintf("HER2 sizes: pos=%d neg=%d", n_pos, n_neg),
  sprintf("Meta template: %s", meta_template_csv),
  sprintf("Meta USED: %s", meta_used_tsv),
  sprintf("Contrast USED: %s", contrast_used_txt),
  sprintf("FDR<0.05: %d", n_fdr),
  sprintf("FDR<0.05 & |log2FC|>=1: %d", n_fc2),
  "Outputs:",
  sprintf("  - %s", file.path(out_dir, "ALL_HER2pos_vs_HER2neg_DE.tsv")),
  sprintf("  - %s", file.path(out_dir, "FILT_HER2pos_vs_HER2neg_DE_FDR0.05.tsv")),
  sprintf("  - %s", file.path(out_dir, "FILT_HER2pos_vs_HER2neg_DE_FDR0.05_LFC1.tsv")),
  sprintf("  - %s", file.path(out_dir, "volcano.png")),
  sprintf("  - %s", file.path(out_dir, "CORR_HER2pos_vs_HER2neg_filtered.png")),
  sprintf("  - %s", file.path(out_dir, "MA.png")),
  pick_line(tt, "ERBB2"),
  pick_line(tt, "GRB7")
)

writeLines(log_txt, file.path(out_dir, "RUN_LOG.txt"))
message("Wolf BP-subtype RNA complete → ", out_dir)
