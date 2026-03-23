# ================= run_brueffer_P76P25_scoring.R =================
# Brueffer (GSE81538) P76/P25 mean-Z scoring using your precomputed overlap lists.
#
# Requires:
#   - Expression matrix: GSE81538_gene_expression_405_transformed.csv
#   - Meta USED: aim2_outputs/GSE81538/GSE81538_RNA_meta_USED.tsv
#   - P76 list CSV: aim2_prioritized/HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv
#   - P25 list CSV: aim2_prioritized/HER2_UP_overlap_25_FDR0.05_LFC1_ALL3_with_stats_Wolf_Robinson_Brueffer.csv
#
# Outputs (in aim2_outputs/GSE81538/):
#   - Brueffer_P76P25_scores_meanZ.tsv
#   - Brueffer_P76_meanZ_boxplot.png
#   - Brueffer_P25_meanZ_boxplot.png
#   - Brueffer_P76_vs_ERBB2_scatter.png (if ERBB2 exists)
#   - Brueffer_RUNLOG_P76P25.txt

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

# ---------------- CONFIG ----------------
base_dir <- "~/Desktop/Thesis/AIM 2 HER2 DATA/brueffer et al"

expr_file <- file.path(base_dir, "GSE81538_gene_expression_405_transformed.csv")
out_dir   <- file.path(base_dir, "aim2_outputs", "GSE81538")
meta_file <- file.path(out_dir, "GSE81538_RNA_meta_USED.tsv")

p76_list_file <- "~/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv"
p25_list_file <- "~/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_25_FDR0.05_LFC1_ALL3_with_stats_Wolf_Robinson_Brueffer.csv"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(expr_file), file.exists(meta_file), file.exists(p76_list_file), file.exists(p25_list_file))

# ---------------- HELPERS ----------------
clean_sym <- function(x) toupper(str_trim(as.character(x)))

z_by_gene <- function(M) {
  mu <- rowMeans(M, na.rm = TRUE)
  sdv <- apply(M, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- NA_real_
  Z <- sweep(M, 1, mu, "-")
  Z <- sweep(Z, 1, sdv, "/")
  Z
}

pick_gene_col <- function(dt) {
  candidates <- c("gene_symbol","Gene","gene","symbol","feature_id","Feature","ID")
  hit <- candidates[candidates %in% names(dt)][1]
  if (!is.na(hit)) return(hit)
  names(dt)[1]
}

# ---------------- READ EXPRESSION ----------------
x <- fread(expr_file, check.names = FALSE)
gcol <- pick_gene_col(x)
setnames(x, gcol, "gene_symbol")

x[, gene_symbol := clean_sym(gene_symbol)]
x <- x[!is.na(gene_symbol) & gene_symbol != ""]
x <- x[!duplicated(gene_symbol)]

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
storage.mode(E) <- "double"

# ---------------- READ META USED ----------------
meta <- fread(meta_file, check.names = FALSE)
stopifnot(all(c("SampleID","Group") %in% names(meta)))
meta$SampleID <- as.character(meta$SampleID)
meta$Group    <- as.character(meta$Group)

keep_samps <- intersect(colnames(E), meta$SampleID)
if (length(keep_samps) < 10) stop("Too few overlapping samples between E and meta USED.")
E <- E[, keep_samps, drop = FALSE]
meta <- meta[match(keep_samps, meta$SampleID)]

# ---------------- READ GENE LISTS ----------------
p76_dt <- fread(p76_list_file, check.names = FALSE)
p25_dt <- fread(p25_list_file, check.names = FALSE)

p76_gene_col <- pick_gene_col(p76_dt)
p25_gene_col <- pick_gene_col(p25_dt)

P76 <- unique(clean_sym(p76_dt[[p76_gene_col]]))
P25 <- unique(clean_sym(p25_dt[[p25_gene_col]]))

P76 <- P76[!is.na(P76) & P76 != ""]
P25 <- P25[!is.na(P25) & P25 != ""]

p76_matched <- intersect(P76, rownames(E))
p25_matched <- intersect(P25, rownames(E))

message("Matched P76 genes: ", length(p76_matched))
message("Matched P25 genes: ", length(p25_matched))

if (length(p76_matched) < 20) stop("Too few P76 genes matched in Brueffer (", length(p76_matched), ").")
if (length(p25_matched) < 10) stop("Too few P25 genes matched in Brueffer (", length(p25_matched), ").")

# ---------------- SCORE (mean Z) ----------------
Z <- z_by_gene(E)

p76_score <- colMeans(Z[p76_matched, , drop = FALSE], na.rm = TRUE)
p25_score <- colMeans(Z[p25_matched, , drop = FALSE], na.rm = TRUE)

out <- data.table(
  SampleID = colnames(E),
  Group    = meta$Group,
  P76_meanZ = as.numeric(p76_score),
  P25_meanZ = as.numeric(p25_score),
  n_P76_genes_matched = length(p76_matched),
  n_P25_genes_matched = length(p25_matched)
)

out_tsv <- file.path(out_dir, "Brueffer_P76P25_scores_meanZ.tsv")
fwrite(out, out_tsv, sep = "\t")
message("Saved: ", out_tsv)

# ---------------- PLOTS ----------------
p_box <- function(df, ycol, ttl, sub) {
  ggplot(df, aes(x = Group, y = .data[[ycol]], fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.45, size = 1.2) +
    theme_classic(base_size = 18) +
    labs(title = ttl, subtitle = sub, x = "Group", y = ycol)
}

g1 <- p_box(out, "P76_meanZ",
            "Brueffer: P76 meanZ by HER2 status",
            paste0("Matched: ", length(p76_matched), "/", length(P76)))
ggsave(file.path(out_dir, "Brueffer_P76_meanZ_boxplot.png"), g1, width = 10, height = 7, dpi = 220)

g2 <- p_box(out, "P25_meanZ",
            "Brueffer: P25 meanZ by HER2 status",
            paste0("Matched: ", length(p25_matched), "/", length(P25)))
ggsave(file.path(out_dir, "Brueffer_P25_meanZ_boxplot.png"), g2, width = 10, height = 7, dpi = 220)

# Optional: correlate P76 with ERBB2 expression
if ("ERBB2" %in% rownames(E)) {
  erbb2 <- as.numeric(E["ERBB2", out$SampleID])
  rho <- suppressWarnings(cor(out$P76_meanZ, erbb2, method = "spearman", use = "complete.obs"))
  g3 <- ggplot(out, aes(x = P76_meanZ, y = erbb2, color = Group)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_classic(base_size = 18) +
    labs(title = "Brueffer: P76 meanZ vs ERBB2 expression",
         subtitle = paste0("Spearman rho = ", signif(rho, 3)),
         x = "P76 meanZ", y = "ERBB2 expression")
  ggsave(file.path(out_dir, "Brueffer_P76_vs_ERBB2_scatter.png"), g3, width = 10, height = 7, dpi = 220)
}

# ---------------- RUN LOG ----------------
log_txt <- c(
  "Brueffer P76/P25 scoring",
  paste0("expr_file: ", expr_file),
  paste0("meta_used: ", meta_file),
  paste0("p76_list: ", p76_list_file),
  paste0("p25_list: ", p25_list_file),
  paste0("samples_scored: ", nrow(out)),
  paste0("P76_matched: ", length(p76_matched)),
  paste0("P25_matched: ", length(p25_matched)),
  paste0("output_tsv: ", out_tsv)
)
writeLines(log_txt, file.path(out_dir, "Brueffer_RUNLOG_P76P25.txt"))

message("DONE: Brueffer P76/P25 scoring complete.")
# ================= end =================
