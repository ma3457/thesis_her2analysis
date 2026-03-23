# ============================
# run_Wolf_scoring.R
#
# This script computes simple gene-set scores for the P76 and P25 gene lists
# in the Wolf (GSE194040 / I-SPY2) RNA expression dataset.
#
# The script:
#  - Loads the Wolf RNA expression matrix
#  - Loads HER2 labels based on BP-subtype classification
#  - Loads HArPS labels derived from RPS-5 annotations
#  - Loads the P76 and P25 gene lists (cross-cohort overlap sets)
#  - Computes per-sample scores using mean z-scored expression
#  - Generates basic plots:
#      * Score vs HER2 status
#      * Score vs HArPS label (within HER2-negative tumors)
#      * P76 vs P25 score correlation
#  - Writes a table of per-sample scores and summary plots to disk
#
# This script does not perform differential expression analysis.
# It only computes and visualizes gene-set-based scores.
# ============================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# ---------------- script directory ----------------
script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) && file.exists(f)) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(".", winslash = "/")
})

# ---------------- paths ----------------
wolf_dir <- script_dir

expr_file <- file.path(
  wolf_dir,
  "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt"
)

bp_lab_file <- file.path(wolf_dir, "labels_from_BP_subtype.csv")

rps5_labels_file <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/labels_from_RPS5.csv"

p76_file <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv"
p25_file <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_25_FDR0.05_LFC1_ALL3_with_stats_Wolf_Robinson_Brueffer.csv"

out_dir <- file.path(wolf_dir, "aim2_outputs", "Wolf_P76P25_scoring")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(expr_file),
  file.exists(bp_lab_file),
  file.exists(rps5_labels_file),
  file.exists(p76_file),
  file.exists(p25_file)
)

# ---------------- helper functions ----------------

wrap_title <- function(s, width = 60) paste(strwrap(s, width = width), collapse = "\n")

read_gene_set <- function(fp) {
  dt <- fread(fp)
  cand <- c("gene_symbol","GeneSymbol","symbol","gene","Gene","feature_id","feature","Feature")
  col  <- cand[cand %in% names(dt)][1]
  if (is.na(col)) stop("No gene column found in ", fp, "\nColumns: ", paste(names(dt), collapse=", "))
  g <- unique(str_trim(as.character(dt[[col]])))
  g[nzchar(g)]
}

row_z <- function(M) {
  mu <- rowMeans(M, na.rm = TRUE)
  sdv <- apply(M, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  (M - mu) / sdv
}

score_mean_z <- function(Ez, gset) {
  if (length(gset) < 2) return(rep(NA_real_, ncol(Ez)))
  colMeans(Ez[gset, , drop = FALSE], na.rm = TRUE)
}

plot_with_margins <- function(expr) {
  op <- par(mar = c(5, 5, 5, 2))
  on.exit(par(op), add = TRUE)
  force(expr)
}

plot_box <- function(y, groups, main, ylab, out_png) {
  png(out_png, width = 1500, height = 1100, res = 160)
  plot_with_margins({
    boxplot(split(y, groups),
            main = wrap_title(main, 65),
            ylab = ylab, xlab = "",
            col = "grey90", border = "grey30")
    stripchart(split(y, groups), vertical = TRUE, method = "jitter",
               pch = 16, cex = 0.6, add = TRUE, col = "grey40")
  })
  dev.off()
}

# ---------------- load gene sets ----------------

P76 <- read_gene_set(p76_file)
P25 <- read_gene_set(p25_file)

# ---------------- read expression matrix ----------------

x <- fread(expr_file, check.names = FALSE)

if (is.na(names(x)[1]) || names(x)[1] == "") setnames(x, 1, "gene_symbol")
if (!names(x)[1] %in% c("gene_symbol", "Gene", "feature_id", "ID", "Feature")) {
  setnames(x, 1, "gene_symbol")
}

x <- x[!is.na(gene_symbol) & str_trim(gene_symbol) != ""]
x <- x[!duplicated(gene_symbol), ]

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
rownames(E) <- make.unique(as.character(x$gene_symbol))
storage.mode(E) <- "numeric"

vals <- as.numeric(E[is.finite(E)])
if (length(vals) && (quantile(vals, 0.99, na.rm = TRUE) > 30 || max(vals, na.rm = TRUE) > 1000)) {
  E <- log2(E + 1)
}

zv <- apply(E, 1, sd, na.rm = TRUE) == 0
if (any(zv)) E <- E[!zv, , drop = FALSE]

# ---------------- read BP-subtype HER2 labels ----------------

bp <- fread(bp_lab_file, check.names = FALSE)
stopifnot(all(c("sample_id", "group") %in% names(bp)))

bp[, sample_id := str_trim(as.character(sample_id))]
bp[, group := str_trim(as.character(group))]

keep <- intersect(colnames(E), bp$sample_id)
if (!length(keep)) stop("No overlap between expression columns and BP-subtype labels.")

E  <- E[, keep, drop = FALSE]
bp <- bp[match(colnames(E), bp$sample_id)]

ok <- !is.na(bp$group) & bp$group %in% c("Positive","Negative")
E  <- E[, ok, drop = FALSE]
bp <- bp[ok]

grp <- factor(ifelse(bp$group == "Positive", "HER2pos", "HER2neg"),
              levels = c("HER2neg", "HER2pos"))

# ---------------- read RPS-5 HArPS labels ----------------

rps <- fread(rps5_labels_file, check.names = FALSE)
stopifnot(all(c("sample_id", "HARP") %in% names(rps)))

rps[, sample_id := str_trim(as.character(sample_id))]
rps[, HARP := str_trim(as.character(HARP))]

harp_vec <- rps$HARP[match(colnames(E), rps$sample_id)]

# ---------------- compute scores ----------------

P76_in <- intersect(P76, rownames(E))
P25_in <- intersect(P25, rownames(E))

Ez <- row_z(E)

P76_meanZ <- score_mean_z(Ez, P76_in)
P25_meanZ <- score_mean_z(Ez, P25_in)

# ---------------- save score table ----------------

score_tbl <- data.table(
  SampleID = colnames(E),
  HER2_group_BP = as.character(grp),
  HARP_RPS5 = harp_vec,
  P76_meanZ = as.numeric(P76_meanZ),
  P25_meanZ = as.numeric(P25_meanZ),
  n_P76_genes_matched = length(P76_in),
  n_P25_genes_matched = length(P25_in)
)

fwrite(score_tbl, file.path(out_dir, "Wolf_P76P25_scores_meanZ.tsv"), sep = "\t")

# ---------------- plots ----------------

plot_box(
  score_tbl$P76_meanZ, grp,
  "Wolf: P76 score by HER2 status (BP-subtype)",
  "P76 score (mean z)",
  file.path(out_dir, "P76_meanZ_by_HER2_BP.png")
)

plot_box(
  score_tbl$P25_meanZ, grp,
  "Wolf: P25 score by HER2 status (BP-subtype)",
  "P25 score (mean z)",
  file.path(out_dir, "P25_meanZ_by_HER2_BP.png")
)

idx <- which(grp == "HER2neg" & !is.na(harp_vec) & harp_vec %in% c("HARP+","HARP-"))
if (length(idx) >= 6) {
  harp_grp <- factor(harp_vec[idx], levels = c("HARP-","HARP+"))
  
  plot_box(
    score_tbl$P76_meanZ[idx], harp_grp,
    "Wolf (HER2- only): P76 score by HArPS label",
    "P76 score (mean z)",
    file.path(out_dir, "P76_meanZ_by_HARP_within_HER2neg.png")
  )
  
  plot_box(
    score_tbl$P25_meanZ[idx], harp_grp,
    "Wolf (HER2- only): P25 score by HArPS label",
    "P25 score (mean z)",
    file.path(out_dir, "P25_meanZ_by_HARP_within_HER2neg.png")
  )
}

png(file.path(out_dir, "P76_vs_P25_meanZ_scatter.png"), width = 1400, height = 1100, res = 160)
plot_with_margins({
  plot(score_tbl$P76_meanZ, score_tbl$P25_meanZ,
       pch = 16, cex = 0.7, col = "grey40",
       xlab = "P76 score (mean z)",
       ylab = "P25 score (mean z)",
       main = wrap_title("Wolf: P76 vs P25 scores", 65))
  abline(lm(score_tbl$P25_meanZ ~ score_tbl$P76_meanZ), lwd = 2)
})
dev.off()

# ---------------- simple log ----------------

log_txt <- c(
  "Wolf P76/P25 scoring run",
  paste0("Expression file: ", expr_file),
  paste0("BP labels: ", bp_lab_file),
  paste0("RPS5 labels: ", rps5_labels_file),
  paste0("P76 genes matched: ", length(P76_in)),
  paste0("P25 genes matched: ", length(P25_in)),
  paste0("Samples scored: ", ncol(E))
)

writeLines(log_txt, file.path(out_dir, "RUN_LOG.txt"))

message("Done. Outputs in: ", out_dir)
