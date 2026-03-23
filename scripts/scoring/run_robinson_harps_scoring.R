# ================= run_robinson_harps_scoring.R =================
# Score Robinson (GSE199633) using a Wolf-derived HARPS signature (TOP N genes by FDR among logFC>0).

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# ---------------- PATHS ----------------
base_dir <- "~/Desktop/Thesis/AIM 2 HER2 DATA/robinson et al"
res_dir  <- file.path(base_dir, "aim2_outputs", "GSE199633_Robinson2025")

expr_file <- file.path(base_dir, "GSE199633_expr_for_pipeline.tsv")
meta_file <- file.path(res_dir,  "GSE199633_Robinson2025_RNA_meta_USED.tsv")

# Wolf HARPS FDR-only table
harps_de_file <- "~/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/Wolf_HARP_RNA/LIMMA_FDR0.05.csv"

stopifnot(dir.exists(base_dir), dir.exists(res_dir))
stopifnot(file.exists(expr_file), file.exists(meta_file))
stopifnot(file.exists(harps_de_file))

# merge with Robinson P76/P25 scores if present
p76p25_file <- file.path(res_dir, "Robinson_P76P25_scores_meanZ.tsv")
has_p76p25  <- file.exists(p76p25_file)

# ---------------- SETTINGS ----------------
TOP_N_HARPS <- 200  # <- defensible HARPS signature size (try 200 or 300; keep fixed across cohorts)

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

# ---------------- READ ROBINSON EXPRESSION ----------------
Xdf <- fread(expr_file, check.names = FALSE)
gene_col <- names(Xdf)[1]

genes <- clean_sym(Xdf[[gene_col]])
Xdf[[gene_col]] <- NULL

X <- as.matrix(Xdf)
rownames(X) <- genes
storage.mode(X) <- "double"

# ---------------- READ META ----------------
meta <- fread(meta_file, check.names = FALSE)
stopifnot(all(c("SampleID","Group") %in% names(meta)))

meta$SampleID <- as.character(meta$SampleID)
meta$Group    <- as.character(meta$Group)

keep_samps <- intersect(colnames(X), meta$SampleID)
if (length(keep_samps) < 10) stop("Too few overlapping samples between expr and meta USED.")

X <- X[, keep_samps, drop = FALSE]
meta <- meta[match(keep_samps, meta$SampleID)]

# ---------------- READ WOLF HARPS DE ----------------
harps_dt <- fread(harps_de_file, check.names = FALSE)

gcol <- detect_gene_col(harps_dt)
fccol <- detect_logfc_col(harps_dt)
qcol <- detect_fdr_col(harps_dt)

harps_dt[, gene := clean_sym(get(gcol))]
harps_dt[, logFC := as.numeric(get(fccol))]
harps_dt[, FDR := as.numeric(get(qcol))]
harps_dt <- harps_dt[!is.na(gene) & gene != "" & is.finite(logFC) & is.finite(FDR)]

message("Using HARPS DE file: ", harps_de_file)
message("Wolf HARPS table rows after cleaning: ", nrow(harps_dt))

# ---------------- BUILD HARPS SIGNATURE (TOP N UP GENES BY FDR) ----------------
harps_up_ranked <- harps_dt[logFC > 0][order(FDR)]
HARPS_UP <- unique(harps_up_ranked[1:min(TOP_N_HARPS, .N), gene])

message("HARPS signature size (TOP_N_HARPS): ", length(HARPS_UP))

# Match to Robinson
harps_matched <- intersect(HARPS_UP, rownames(X))
message("HARPS genes matched in Robinson: ", length(harps_matched))

if (length(harps_matched) < 20) {
  stop("Too few HARPS genes matched in Robinson (", length(harps_matched),
       "). Check gene symbol mapping / expression rownames.")
}

# ---------------- SCORE (mean Z) ----------------
Z <- z_by_gene(X)
harps_score <- colMeans(Z[harps_matched, , drop = FALSE], na.rm = TRUE)

out <- data.table(
  SampleID = names(harps_score),
  Group    = meta$Group,
  HARPS_meanZ = as.numeric(harps_score),
  n_HARPS_genes_matched = length(harps_matched),
  TOP_N_HARPS = TOP_N_HARPS
)

# ---------------- SAVE ----------------
out_tsv <- file.path(res_dir, "Robinson_HARPS_scores_meanZ.tsv")
fwrite(out, out_tsv, sep = "\t")
message("Saved: ", out_tsv)

# ---------------- PLOTS ----------------
p1 <- ggplot(out, aes(x = Group, y = HARPS_meanZ, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.2) +
  labs(
    title = "Robinson: HARPS meanZ by HER2 status",
    subtitle = paste0("TOP_N_HARPS=", TOP_N_HARPS, " | matched=", length(harps_matched)),
    x = "Group", y = "HARPS meanZ"
  ) +
  theme_classic(base_size = 18)

ggsave(file.path(res_dir, "Robinson_HARPS_meanZ_boxplot_by_HER2.png"),
       p1, width = 10, height = 7, dpi = 200)

out_neg <- out %>% filter(Group == "HER2neg")
p2 <- ggplot(out_neg, aes(x = HARPS_meanZ)) +
  geom_histogram(bins = 35) +
  labs(
    title = "Robinson: HARPS meanZ distribution within HER2− tumors",
    subtitle = paste0("n(HER2−) = ", nrow(out_neg)),
    x = "HARPS meanZ", y = "Count"
  ) +
  theme_classic(base_size = 18)

ggsave(file.path(res_dir, "Robinson_HARPS_meanZ_hist_HER2neg.png"),
       p2, width = 10, height = 6, dpi = 200)

# ---------------- MERGE WITH P76/P25 ----------------
if (has_p76p25) {
  p76p25 <- fread(p76p25_file, check.names = FALSE)
  p76p25$SampleID <- as.character(p76p25$SampleID)
  
  joint <- left_join(p76p25, out, by = "SampleID")
  joint_out <- file.path(res_dir, "Robinson_joint_scores_P76P25_HARPS.tsv")
  fwrite(as.data.table(joint), joint_out, sep = "\t")
  message("Saved: ", joint_out)
  
  if ("P76_meanZ" %in% names(joint)) {
    p3 <- ggplot(joint, aes(x = P76_meanZ, y = HARPS_meanZ, color = Group)) +
      geom_point(alpha = 0.7, size = 2) +
      labs(
        title = "Robinson: P76 vs HARPS (two-axis HER2 framework)",
        subtitle = paste0("TOP_N_HARPS=", TOP_N_HARPS),
        x = "P76 meanZ (amplification axis)",
        y = "HARPS meanZ (HER2-like axis)"
      ) +
      theme_classic(base_size = 18)
    
    ggsave(file.path(res_dir, "Robinson_joint_P76_vs_HARPS.png"),
           p3, width = 10, height = 7, dpi = 200)
  }
}

message("DONE: Robinson HARPS scoring complete.")
# ================= end =================
