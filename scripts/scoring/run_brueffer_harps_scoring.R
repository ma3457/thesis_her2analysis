# ================= run_brueffer_harps_scoring.R =================
# Brueffer (GSE81538) HARPS mean-Z scoring using Wolf-derived HARPS signature (TOP N genes by FDR among logFC>0).
#
# Requires:
#   - Expression: GSE81538_gene_expression_405_transformed.csv
#   - Meta USED: aim2_outputs/GSE81538/GSE81538_RNA_meta_USED.tsv
#   - Wolf HARPS DE: wolf et al/aim2_outputs/Wolf_HARP_RNA/LIMMA_FDR0.05.csv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# ---------------- CONFIG ----------------
base_dir <- "~/Desktop/Thesis/AIM 2 HER2 DATA/brueffer et al"

expr_file <- file.path(base_dir, "GSE81538_gene_expression_405_transformed.csv")
out_dir   <- file.path(base_dir, "aim2_outputs", "GSE81538")
meta_file <- file.path(out_dir, "GSE81538_RNA_meta_USED.tsv")

harps_de_file <- "~/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/Wolf_HARP_RNA/LIMMA_FDR0.05.csv"
TOP_N_HARPS   <- 200

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(expr_file), file.exists(meta_file), file.exists(harps_de_file))

# ---------------- HELPERS ----------------
clean_sym <- function(x) toupper(trimws(as.character(x)))

z_by_gene <- function(M) {
  mu <- rowMeans(M, na.rm = TRUE)
  sdv <- apply(M, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- NA_real_
  Z <- sweep(M, 1, mu, "-")
  Z <- sweep(Z, 1, sdv, "/")
  Z
}

detect_gene_col <- function(dt) {
  candidates <- c("gene_symbol","Gene","gene","symbol","feature_id","ID","id")
  hit <- candidates[candidates %in% names(dt)][1]
  if (!is.na(hit)) return(hit)
  names(dt)[1]
}
detect_logfc_col <- function(dt) {
  candidates <- c("logFC","log2FC","log2FoldChange","log2fc")
  hit <- candidates[candidates %in% names(dt)][1]
  if (!is.na(hit)) return(hit)
  hit2 <- names(dt)[grepl("log", names(dt), ignore.case = TRUE) &
                      grepl("fc",  names(dt), ignore.case = TRUE)][1]
  if (is.na(hit2)) stop("Could not find logFC column in HARPS DE file.")
  hit2
}
detect_fdr_col <- function(dt) {
  candidates <- c("adj.P.Val","FDR","padj","qvalue")
  hit <- candidates[candidates %in% names(dt)][1]
  if (!is.na(hit)) return(hit)
  hit2 <- names(dt)[grepl("adj", names(dt), ignore.case = TRUE) &
                      grepl("p",   names(dt), ignore.case = TRUE)][1]
  if (is.na(hit2)) stop("Could not find FDR/adj.P.Val column in HARPS DE file.")
  hit2
}

pick_gene_col <- function(dt) {
  candidates <- c("gene_symbol","Gene","feature_id","ID","Feature")
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

# ---------------- META USED ----------------
meta <- fread(meta_file, check.names = FALSE)
stopifnot(all(c("SampleID","Group") %in% names(meta)))
meta$SampleID <- as.character(meta$SampleID)
meta$Group    <- as.character(meta$Group)

keep_samps <- intersect(colnames(E), meta$SampleID)
if (length(keep_samps) < 10) stop("Too few overlapping samples between E and meta USED.")
E <- E[, keep_samps, drop = FALSE]
meta <- meta[match(keep_samps, meta$SampleID)]

# ---------------- READ WOLF HARPS DE ----------------
harps_dt <- fread(harps_de_file, check.names = FALSE)

gcol2 <- detect_gene_col(harps_dt)
fccol <- detect_logfc_col(harps_dt)
qcol  <- detect_fdr_col(harps_dt)

harps_dt[, gene := clean_sym(get(gcol2))]
harps_dt[, logFC := as.numeric(get(fccol))]
harps_dt[, FDR := as.numeric(get(qcol))]
harps_dt <- harps_dt[!is.na(gene) & gene != "" & is.finite(logFC) & is.finite(FDR)]

message("Using HARPS DE file: ", harps_de_file)
message("Wolf HARPS table rows after cleaning: ", nrow(harps_dt))

harps_up_ranked <- harps_dt[logFC > 0][order(FDR)]
HARPS_UP <- unique(harps_up_ranked[1:min(TOP_N_HARPS, .N), gene])

message("HARPS signature size (TOP_N_HARPS): ", length(HARPS_UP))

harps_matched <- intersect(HARPS_UP, rownames(E))
message("HARPS genes matched in Brueffer: ", length(harps_matched))

if (length(harps_matched) < 20) stop("Too few HARPS genes matched in Brueffer (", length(harps_matched), ").")

# ---------------- SCORE ----------------
Z <- z_by_gene(E)
harps_score <- colMeans(Z[harps_matched, , drop = FALSE], na.rm = TRUE)

out <- data.table(
  SampleID = colnames(E),
  Group    = meta$Group,
  HARPS_meanZ = as.numeric(harps_score),
  n_HARPS_genes_matched = length(harps_matched),
  TOP_N_HARPS = TOP_N_HARPS
)

out_tsv <- file.path(out_dir, "Brueffer_HARPS_scores_meanZ.tsv")
fwrite(out, out_tsv, sep = "\t")
message("Saved: ", out_tsv)

# ---------------- PLOTS ----------------
p1 <- ggplot(out, aes(x = Group, y = HARPS_meanZ, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.2) +
  labs(
    title = "Brueffer: HARPS meanZ by HER2 status",
    subtitle = paste0("TOP_N_HARPS=", TOP_N_HARPS, " | matched=", length(harps_matched)),
    x = "Group", y = "HARPS meanZ"
  ) +
  theme_classic(base_size = 18)

ggsave(file.path(out_dir, "Brueffer_HARPS_meanZ_boxplot_by_HER2.png"),
       p1, width = 10, height = 7, dpi = 220)

out_neg <- out %>% filter(Group == "HER2neg")
p2 <- ggplot(out_neg, aes(x = HARPS_meanZ)) +
  geom_histogram(bins = 35) +
  labs(
    title = "Brueffer: HARPS meanZ distribution within HER2− tumors",
    subtitle = paste0("n(HER2−) = ", nrow(out_neg)),
    x = "HARPS meanZ", y = "Count"
  ) +
  theme_classic(base_size = 18)

ggsave(file.path(out_dir, "Brueffer_HARPS_meanZ_hist_HER2neg.png"),
       p2, width = 10, height = 6, dpi = 220)

# ---------------- RUN LOG ----------------
log_txt <- c(
  "Brueffer HARPS scoring",
  paste0("expr_file: ", expr_file),
  paste0("meta_used: ", meta_file),
  paste0("harps_de_file: ", harps_de_file),
  paste0("TOP_N_HARPS: ", TOP_N_HARPS),
  paste0("samples_scored: ", nrow(out)),
  paste0("HARPS_matched: ", length(harps_matched)),
  paste0("output_tsv: ", out_tsv)
)
writeLines(log_txt, file.path(out_dir, "Brueffer_RUNLOG_HARPS.txt"))

message("DONE: Brueffer HARPS scoring complete.")
# ================= end =================
