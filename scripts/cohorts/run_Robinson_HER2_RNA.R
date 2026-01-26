# ================= run_Robinson_HER2_RNA.R — Robinson RNA (GSE199633) =================
# Purpose:
#   - Build a gene-level matrix from the GEO series_matrix + GPL mapping
#   - Run her2_pipeline.R (limma DE) using study-defined HER2 status
#   - Write meta artifacts (template + used + contrast) for reproducibility
#
# Outputs (in aim2_outputs/GSE199633_Robinson2025):
#   - *_DE_HER2pos_vs_HER2neg.tsv              (FDR<0.05 hits from her2_pipeline)
#   - FILT_*_FDR0.05_LFC1.tsv                  (subset: |log2FC|>=1)
#   - CORR_*_filtered.png (or fallback)        (sample correlation heatmap)
#   - *_RNA_meta_template.csv                  (editable metadata snapshot)
#   - *_RNA_meta_USED.tsv                      (samples actually modeled)
#   - *_RNA_contrast_USED.txt                  (contrast used)

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(GEOquery)
  library(dplyr)
  library(pheatmap)
})

# ---- self-locating (no setwd required) ----
script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) && file.exists(f)) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(".", winslash = "/")
})

# ---------------- CONFIG ----------------
base_dir <- script_dir
series   <- file.path(base_dir, "GSE199633_series_matrix.txt")
stopifnot(file.exists(series))

# her2_pipeline.R is expected to live one level up from cohort folders
pipe_path <- file.path(dirname(base_dir), "her2_pipeline.R")
stopifnot(file.exists(pipe_path))

out_prefix   <- "GSE199633_Robinson2025"
series_files <- c("GSE199633_series_matrix.txt")

res_dir <- file.path(base_dir, "aim2_outputs", out_prefix)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# meta artifacts
meta_template_csv <- file.path(base_dir, paste0(out_prefix, "_RNA_meta_template.csv"))
meta_used_tsv     <- file.path(res_dir,  paste0(out_prefix, "_RNA_meta_USED.tsv"))
contrast_used_txt <- file.path(res_dir,  paste0(out_prefix, "_RNA_contrast_USED.txt"))

# helper mapping file
gsm_to_title_tsv  <- file.path(base_dir, paste0(out_prefix, "_GSM_to_title.tsv"))

# ---------------- 1) Extract expression block from SERIES ----------------
L  <- readLines(series, warn = FALSE)
i1 <- which(L == "!series_matrix_table_begin")
i2 <- which(L == "!series_matrix_table_end")
if (!length(i1) || !length(i2) || i2 <= i1) stop("No matrix block in series_matrix file.")

tmp <- tempfile(fileext = ".txt")
writeLines(L[(i1 + 1):(i2 - 1)], tmp)
tab <- fread(tmp, data.table = FALSE, check.names = FALSE)
unlink(tmp)

grab <- function(rx) {
  ln <- grep(rx, L, value = TRUE)[1]
  if (length(ln) == 0 || is.na(ln)) stop("Missing field: ", rx)
  gsub('^"|"$', "", strsplit(ln, "\t", fixed = TRUE)[[1]][-1])
}

gsm   <- grab("^!Sample_geo_accession")
title <- grab("^!Sample_title")
stopifnot(length(gsm) == length(title))

# GSM->title mapping (handy for meta exports)
gmap <- data.table(GEO_accession = gsm, Sample_title = title)
fwrite(gmap, gsm_to_title_tsv, sep = "\t", quote = FALSE, na = "NA")
message("Wrote GSM->title map: ", gsm_to_title_tsv)

pick_feat <- function(nm) {
  nm_low <- tolower(nm)
  prefs <- c("id_ref","id","genesymbol","gene_symbol","symbol","gene","ensembl","ensembl_gene_id")
  for (p in prefs) {
    hit <- which(nm_low == p)
    if (length(hit)) return(hit[1])
  }
  1L
}

fcol <- pick_feat(colnames(tab))
colnames(tab)[fcol] <- "feature_id"

# matrix block sometimes lacks GSM headers; fallback to positional selection if needed
gsm_cols <- intersect(colnames(tab), gsm)
if (length(gsm_cols) < 6) {
  message("WARNING: GSM headers not found in matrix block; using last N columns by position.")
  if (ncol(tab) < (1 + length(gsm))) stop("Matrix block does not have enough columns for samples.")
  gsm_cols <- tail(colnames(tab), length(gsm))
}

expr <- cbind(feature_id = tab[["feature_id"]], tab[, gsm_cols, drop = FALSE])

# numeric coercion
numcols <- setdiff(colnames(expr), "feature_id")
for (cn in numcols) {
  v <- expr[[cn]]
  if (!is.numeric(v)) v <- suppressWarnings(as.numeric(as.character(v)))
  expr[[cn]] <- v
}

keep <- rowSums(is.finite(as.matrix(expr[, numcols, drop = FALSE]))) > 0
expr <- expr[keep, , drop = FALSE]
message("Probe-level matrix: rows=", nrow(expr), "  samples=", length(numcols))

# ---------------- 2) Probe -> gene_symbol via dominant GPL ----------------
splat <- grab("^!Sample_platform_id")
major_gpl <- names(which.max(table(splat)))
message("Using platform: ", major_gpl)

gpl <- getGEO(major_gpl)
gt  <- Table(gpl)

sym_col <- names(gt)[tolower(names(gt)) %in% c(
  "gene symbol","gene_symbol","symbol","genesymbol","hugo symbol","hgnc symbol","hgnc_symbol"
)]
if (!length(sym_col) && "Gene Assignment" %in% names(gt)) sym_col <- "Gene Assignment"
if (!length(sym_col)) stop("No gene symbol column detected in GPL table.")
sym_col <- sym_col[1]

map <- gt[, c("ID", sym_col)]
colnames(map) <- c("probe_id","sym_raw")

norm_sym <- function(x) {
  x <- as.character(x)
  x <- gsub(" ///.*$", "", x)
  x <- gsub("\\s*\\(.*\\)$", "", x)
  x <- gsub("^LOC[0-9]+$", "", x)
  x <- trimws(x)
  x[nchar(x) == 0] <- NA
  x
}

map$gene_symbol <- norm_sym(map$sym_raw)
map <- map[!is.na(map$gene_symbol) & map$gene_symbol != "", c("probe_id","gene_symbol")]

expr2 <- dplyr::inner_join(expr, map, by = c("feature_id" = "probe_id"))
if (nrow(expr2) < 2000) stop("Too few mapped probes after GPL join: ", nrow(expr2))

gene_expr <- expr2 %>%
  dplyr::select(gene_symbol, dplyr::all_of(numcols)) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarize(dplyr::across(dplyr::everything(), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

# force sample columns to GSM IDs
sample_cols <- setdiff(colnames(gene_expr), "gene_symbol")
if (length(sample_cols) != length(gsm)) {
  stop("Sample column count mismatch: gene_expr has ", length(sample_cols),
       " sample cols but GSM list has ", length(gsm),
       ". This series_matrix needs manual inspection.")
}
setnames(gene_expr, old = sample_cols, new = gsm)
message("PRE-FLIGHT OK: gene_expr columns are GSM IDs.")

expr_gene_tsv <- file.path(base_dir, "GSE199633_expr_for_pipeline.tsv")
fwrite(gene_expr, expr_gene_tsv, sep = "\t", quote = FALSE, na = "NA")
message("Gene-level matrix written: ", expr_gene_tsv, "  genes=", nrow(gene_expr))

# ---------------- 3) Run her2_pipeline ----------------
source(pipe_path)
options(her2.transform_mode = "auto")

obj <- her2_run(
  dataset_dir         = base_dir,
  expr_file           = basename(expr_gene_tsv),
  series_files        = series_files,
  out_prefix          = out_prefix,
  collapse_replicates = FALSE,
  pos_label           = "Positive",
  neg_label           = "Negative",
  her2_definition     = "study-defined"
)

write_run_log(base_dir, out_prefix)

# ---------------- 3b) META TEMPLATE ----------------
if (!file.exists(meta_template_csv)) {
  g_source <- tryCatch(grab("^!Sample_source_name_ch1"), error = function(e) rep(NA_character_, length(gsm)))
  g_ch1    <- tryCatch(grab("^!Sample_characteristics_ch1"), error = function(e) rep(NA_character_, length(gsm)))
  
  meta_template <- data.table(
    GEO_accession = gsm,
    Sample_title  = title,
    Source_name   = g_source,
    Characteristics_ch1 = g_ch1,
    Notes = ""
  )
  fwrite(meta_template, meta_template_csv)
  message("Wrote meta template: ", meta_template_csv)
} else {
  message("Meta template exists: ", meta_template_csv)
}

# ---------------- 4) META USED + CONTRAST USED ----------------
stopifnot(!is.null(obj$E), !is.null(obj$ann))

used_ids <- colnames(obj$E)
ann_status <- as.character(obj$ann$HER2_status)

Group <- ifelse(ann_status == "Positive", "HER2pos",
                ifelse(ann_status == "Negative", "HER2neg", NA_character_))

meta_used_out <- data.table(
  SampleID = used_ids,
  Group    = Group
)

# attach title if available
if (file.exists(gsm_to_title_tsv)) {
  gmap2 <- fread(file = gsm_to_title_tsv, sep = "\t", header = TRUE)
  if (!"GEO_accession" %in% names(gmap2)) setnames(gmap2, names(gmap2)[1], "GEO_accession")
  if (!"Sample_title"  %in% names(gmap2)) setnames(gmap2, names(gmap2)[2], "Sample_title")
  
  meta_used_out[, GEO_accession := SampleID]
  meta_used_out <- merge(meta_used_out, gmap2, by = "GEO_accession", all.x = TRUE)
  setcolorder(meta_used_out, c("SampleID","Group","GEO_accession","Sample_title"))
}

fwrite(meta_used_out, meta_used_tsv, sep = "\t", quote = FALSE, na = "NA")
message("Wrote meta USED: ", meta_used_tsv)

n_pos <- sum(meta_used_out$Group == "HER2pos", na.rm = TRUE)
n_neg <- sum(meta_used_out$Group == "HER2neg", na.rm = TRUE)

writeLines(c(
  "Mode: HER2",
  "Contrast: HER2pos - HER2neg",
  paste0("Cohort: ", out_prefix),
  paste0("n_pos: ", n_pos),
  paste0("n_neg: ", n_neg)
), contrast_used_txt)
message("Wrote contrast USED: ", contrast_used_txt)

# ---------------- 5) FC2 subset + correlation heatmap ----------------
# Note: pipeline DE output is already FDR<0.05 (by design in her2_pipeline.R).
f_in <- file.path(res_dir, paste0(out_prefix, "_DE_HER2pos_vs_HER2neg.tsv"))
stopifnot(file.exists(f_in))

DT <- fread(file = f_in, check.names = FALSE)
stopifnot(all(c("gene_symbol","P.Value","adj.P.Val") %in% names(DT)))

fc_col <- if ("log2FC_HER2pos_vs_HER2neg" %in% names(DT)) "log2FC_HER2pos_vs_HER2neg" else "logFC"
DT[, logFC := as.numeric(get(fc_col))]
DT[, gene_symbol := toupper(str_trim(as.character(gene_symbol)))]

FILT <- DT[abs(logFC) >= 1.0]
fdr_fc2_out <- file.path(res_dir, "FILT_HER2pos_vs_HER2neg_DE_FDR0.05_LFC1.tsv")
fwrite(FILT[, .(gene_symbol, logFC, P.Value, adj.P.Val)], fdr_fc2_out, sep = "\t")
message("Wrote FC2 subset: ", fdr_fc2_out)

# matrix for correlation heatmap
Xdf <- fread(file = expr_gene_tsv)
rn <- Xdf[[1]]; Xdf[[1]] <- NULL
rownames(Xdf) <- rn
X <- as.matrix(Xdf)
storage.mode(X) <- "double"

make_heat <- function(genes, file, ttl) {
  keepg <- intersect(genes, rownames(X))
  if (length(keepg) >= 3) {
    Msub <- X[keepg, , drop = FALSE]
    cm <- cor(t(Msub), use = "pairwise.complete.obs")
    png(file.path(res_dir, file), width = 1400, height = 1200, res = 180)
    pheatmap(cm, main = ttl)
    dev.off()
    TRUE
  } else FALSE
}

heat_ok <- (nrow(FILT) >= 3) && make_heat(
  FILT[["gene_symbol"]],
  "CORR_HER2pos_vs_HER2neg_filtered.png",
  "Correlation — HER2pos vs HER2neg (FDR<0.05 & |log2FC|>=1)"
)

if (!heat_ok) {
  cand <- DT[order(adj.P.Val)][seq_len(min(25, nrow(DT)))]
  make_heat(cand[["gene_symbol"]], "CORR_FALLBACK_top_by_FDR.png",
            "Correlation — Fallback (top by FDR)")
}

message("DONE: ", res_dir)
# ================= end =================
