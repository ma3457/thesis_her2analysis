# ================= run_mertins_PROT.R =================
# GLOBAL Proteome DEP with LIMMA (HER2+ vs HER2−) — Mertins CPTAC BRCA
# Matches RajKumar/Krug conventions:
#   - outputs_root + out_dir_tag(tag)
#   - meta template written once at outputs_root
#   - meta USED + contrast USED (per tag)
#   - DE_full.tsv + filtered (FDR / FDR+FC2) + canonical copies
#   - volcano + MA + QC plots (ERBB2/GRB7) + RUN_LOG
#   - volcano uses -log10(FDR) (adj.P.Val)
# ======================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(limma)
  library(stringr)
})

# ---------------- CONFIG ----------------
base_dir <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/mertins et al"
supp_dir <- file.path(base_dir, "nature18003-s2")

supp1_xlsx <- file.path(supp_dir, "CPTAC_BC_SupplementaryTable01.xlsx")  # sample annotation
supp3_xlsx <- file.path(supp_dir, "CPTAC_BC_SupplementaryTable03.xlsx")  # global proteome matrix

stopifnot(file.exists(supp1_xlsx), file.exists(supp3_xlsx))

# Krug-style naming
out_prefix   <- "Mertins_CPTAC_BRCA_PROT_GLOBAL"
outputs_root <- file.path(base_dir, "aim2_outputs")
dir.create(outputs_root, showWarnings = FALSE)

# Single definition tag for this cohort (clinical HER2 status from Supp01)
TAG <- "ClinicalOnly"
out_dir_tag <- function(tag) file.path(outputs_root, paste0(out_prefix, "_", tag))

# thresholds
min_per_group <- 3
USE_QC_PASS   <- TRUE
DROP_EQUIV    <- TRUE
FC_MIN_LOG2   <- 1
FDR_MAX       <- 0.05

`%||%` <- function(a,b) if (!is.null(a)) a else b

# ---------------- helpers ----------------
extract_tcga <- function(x){
  xU <- toupper(as.character(x))
  out <- rep(NA_character_, length(xU))
  
  m1 <- str_extract(xU, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
  out[!is.na(m1)] <- m1[!is.na(m1)]
  
  need <- is.na(out)
  if (any(need)) {
    xx <- xU[need]
    parts <- str_split(xx, "[\\.-]")
    p1 <- vapply(parts, function(v) if (length(v) >= 2) v[1] else NA_character_, character(1))
    p2 <- vapply(parts, function(v) if (length(v) >= 2) v[2] else NA_character_, character(1))
    ends_tcga <- grepl("TCGA$", xx, ignore.case = TRUE)
    out[need] <- ifelse(ends_tcga & !is.na(p1) & !is.na(p2),
                        paste0("TCGA-", p1, "-", p2),
                        NA_character_)
  }
  out
}

read_matrix_xlsx <- function(path){
  x <- as.data.frame(read_xlsx(path, sheet = 1))
  colnames(x) <- make.unique(colnames(x))
  
  raw_cols <- colnames(x)
  samp_idx <- which(grepl("TCGA", raw_cols, ignore.case = TRUE))
  if (!length(samp_idx)) stop("Could not detect TCGA sample columns in: ", basename(path))
  
  ann_cols  <- raw_cols[-samp_idx]
  samp_cols <- raw_cols[samp_idx]
  
  pick_first <- function(cands) {
    hit <- intersect(cands, ann_cols)
    if (length(hit)) hit[1] else NA_character_
  }
  
  id_primary  <- pick_first(c("geneName","Gene","Gene.Symbol","GeneSymbol","Symbol"))
  id_fallback <- pick_first(c("proteinName","Protein","Name","Description","accession","Accession"))
  
  if (is.na(id_primary) && is.na(id_fallback))
    stop("No usable feature identifier column found in: ", basename(path))
  
  feat <- rep("", nrow(x))
  
  if (!is.na(id_primary)) {
    tmp <- as.character(x[[id_primary]])
    tmp[is.na(tmp)] <- ""
    feat <- tmp
  }
  
  if (!is.na(id_fallback)) {
    fb <- as.character(x[[id_fallback]])
    fb[is.na(fb)] <- ""
    need <- trimws(feat) == ""
    feat[need] <- fb[need]
  }
  
  feat <- trimws(toupper(feat))
  feat <- sub("[,;/].*$", "", feat)          # keep first token if multiple
  feat <- gsub("[^A-Z0-9_.-]", "", feat)     # clean
  feat[feat %in% c("", "NA")] <- ""          # treat literal "NA" as missing
  
  bad <- feat == ""
  if (any(bad)) feat[bad] <- paste0("UNMAPPED_", seq_len(sum(bad)))
  feat <- make.unique(feat)
  
  M <- as.matrix(x[, samp_cols, drop = FALSE])
  storage.mode(M) <- "numeric"
  rownames(M) <- feat
  
  list(M = M, feature_col = paste(c(id_primary, id_fallback), collapse = " | "), sample_cols = samp_cols)
}

robust_impute_and_filter <- function(M){
  all_bad <- rowSums(!is.finite(M)) == ncol(M)
  if (any(all_bad)) M <- M[!all_bad, , drop = FALSE]
  
  M[!is.finite(M)] <- NA
  M[M < -1e6] <- NA
  
  if (anyNA(M)) {
    meds <- apply(M, 1, function(v){
      m <- suppressWarnings(median(v, na.rm = TRUE))
      if (!is.finite(m)) 0 else m
    })
    for (i in seq_len(nrow(M))) if (anyNA(M[i, ])) M[i, is.na(M[i, ])] <- meds[i]
  }
  
  M[!is.finite(M)] <- 0
  sdv <- apply(M, 1, function(v){
    s <- suppressWarnings(sd(v, na.rm = TRUE))
    if (!is.finite(s)) 0 else s
  })
  M <- M[sdv > 0, , drop = FALSE]
  M
}

collapse_duplicate_tcga_cols <- function(M, tcga_ids){
  u <- unique(tcga_ids)
  if (length(u) == length(tcga_ids)) return(list(M = M, tcga = tcga_ids, collapsed = FALSE))
  
  M2 <- sapply(u, function(id){
    cols <- which(tcga_ids == id)
    if (length(cols) == 1) return(M[, cols])
    rowMeans(M[, cols, drop = FALSE], na.rm = TRUE)
  })
  M2 <- as.matrix(M2); rownames(M2) <- rownames(M)
  list(M = M2, tcga = u, collapsed = TRUE)
}

qc_box <- function(E, g_factor, tag, gene){
  gene <- toupper(gene)
  if (!(gene %in% rownames(E))) return(invisible(FALSE))
  out_dir <- out_dir_tag(tag); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fn_pdf <- file.path(out_dir, paste0("QC_", gene, "_boxplot.pdf"))
  
  vals <- as.numeric(E[gene, ])
  v_pos <- vals[g_factor == "HER2pos"]
  v_neg <- vals[g_factor == "HER2neg"]
  m_pos <- mean(v_pos, na.rm = TRUE); s_pos <- sd(v_pos, na.rm = TRUE)
  m_neg <- mean(v_neg, na.rm = TRUE); s_neg <- sd(v_neg, na.rm = TRUE)
  log2fc <- m_pos - m_neg
  ttst <- try(stats::t.test(v_pos, v_neg, var.equal = FALSE), silent = TRUE)
  pval <- if (inherits(ttst, "try-error")) NA_real_ else ttst$p.value
  
  pdf(fn_pdf, 7.2, 6)
  par(mar = c(5,5,4.5,2))
  boxplot(split(vals, g_factor),
          main = sprintf("%s QC — %s", gene, tag),
          xlab = "HER2 group", ylab = paste0(gene, " (expression)"),
          col = "grey90", border = "grey30")
  stripchart(split(vals, g_factor), vertical = TRUE, method = "jitter",
             add = TRUE, pch = 21, bg = "white", col = "black")
  grid(col = "grey90")
  mtext(sprintf("n(pos)=%d n(neg)=%d  mean±SD(pos)=%.2f±%.2f  mean±SD(neg)=%.2f±%.2f  log2FC=%.2f  t-test P=%.3g",
                length(v_pos), length(v_neg), m_pos, s_pos, m_neg, s_neg, log2fc, pval),
        side = 3, line = 0.6, cex = 0.85)
  dev.off()
  TRUE
}

# ---------------- META TEMPLATE (written once) ----------------
meta_template_csv <- file.path(outputs_root, paste0(out_prefix, "_meta_template.csv"))

# ---------------- 1) LOAD proteome matrix (Supp 03) ----------------
mx <- read_matrix_xlsx(supp3_xlsx)
M  <- mx$M

tcga <- extract_tcga(colnames(M))
if (any(is.na(tcga))) {
  bad <- which(is.na(tcga))[1:min(15, sum(is.na(tcga)))]
  stop("Could not extract TCGA IDs from some matrix column names. Example bad cols:\n",
       paste(colnames(M)[bad], collapse = "\n"))
}

cc <- collapse_duplicate_tcga_cols(M, tcga)
M   <- cc$M
colnames(M) <- cc$tcga

drop_unmapped <- grepl("^UNMAPPED_", rownames(M))
if (any(drop_unmapped)) M <- M[!drop_unmapped, , drop = FALSE]

# ---------------- 2) LOAD sample annotation (Supp 01) ----------------
ann_raw <- as.data.frame(read_xlsx(supp1_xlsx, sheet = 1))
colnames(ann_raw) <- make.unique(make.names(colnames(ann_raw)))

tcga_col <- names(ann_raw)[grepl("TCGA", names(ann_raw), ignore.case = TRUE)][1]
her2_col <- names(ann_raw)[grepl("^HER2", names(ann_raw), ignore.case = TRUE) &
                             grepl("Status", names(ann_raw), ignore.case = TRUE)][1]
qc_col   <- names(ann_raw)[grepl("QC", names(ann_raw), ignore.case = TRUE) &
                             grepl("Status", names(ann_raw), ignore.case = TRUE)][1]

if (any(is.na(c(tcga_col, her2_col)))) {
  stop("Could not find TCGA and/or HER2 Status columns in Supp01. Detected:\n",
       "tcga_col=", tcga_col %||% "NA", "  her2_col=", her2_col %||% "NA",
       "  qc_col=", qc_col %||% "NA")
}

clin <- data.frame(
  SampleID = toupper(trimws(as.character(ann_raw[[tcga_col]]))),
  HER2_raw = trimws(as.character(ann_raw[[her2_col]])),
  QC_raw   = if (!is.na(qc_col)) trimws(as.character(ann_raw[[qc_col]])) else NA_character_,
  stringsAsFactors = FALSE
)

clin$Group <- NA_character_
clin$Group[grepl("^pos",  clin$HER2_raw, ignore.case = TRUE)] <- "HER2pos"
clin$Group[grepl("^neg",  clin$HER2_raw, ignore.case = TRUE)] <- "HER2neg"
clin$Group[grepl("equiv", clin$HER2_raw, ignore.case = TRUE)] <- "HER2equiv"

if (USE_QC_PASS && !all(is.na(clin$QC_raw))) {
  clin <- clin[grepl("pass", clin$QC_raw, ignore.case = TRUE), , drop = FALSE]
}
if (DROP_EQUIV) {
  clin <- clin[clin$Group %in% c("HER2pos","HER2neg"), , drop = FALSE]
}

# ---------------- 3) ALIGN ----------------
keep <- intersect(colnames(M), clin$SampleID)
if (!length(keep)) stop("0 overlap between proteome columns and Supp01 TCGA IDs after QC/HER2 filtering.")

M <- M[, keep, drop = FALSE]
clin <- clin[match(colnames(M), clin$SampleID), , drop = FALSE]
stopifnot(all(colnames(M) == clin$SampleID))

grp <- factor(clin$Group, levels = c("HER2neg","HER2pos"))
if (any(is.na(grp))) stop("Missing/invalid HER2 group after alignment.")
if (min(table(grp)) < min_per_group) stop("HER2 groups too small: ", paste(table(grp), collapse=" / "))

# ---------------- 4) META TEMPLATE (create once; allow edits) ----------------
if (!file.exists(meta_template_csv)) {
  meta_tmpl <- data.table(
    SampleID = clin$SampleID,
    HER2_raw = clin$HER2_raw,
    QC_raw   = clin$QC_raw,
    Group_ClinicalOnly = as.character(grp),
    Notes = ""
  )
  fwrite(meta_tmpl, meta_template_csv)
}

# Re-read template to allow manual overrides
meta <- fread(meta_template_csv, data.table = FALSE)
meta$SampleID <- toupper(trimws(as.character(meta$SampleID)))

gcol <- if ("Group_ClinicalOnly" %in% names(meta)) "Group_ClinicalOnly" else "Group"
if (!(gcol %in% names(meta))) stop("Meta template missing Group column. Expected Group_ClinicalOnly.")

tmp <- toupper(trimws(as.character(meta[[gcol]])))
tmp[tmp %in% c("HER2POS","POS","P")] <- "HER2pos"
tmp[tmp %in% c("HER2NEG","NEG","N")] <- "HER2neg"
tmp[!tmp %in% c("HER2pos","HER2neg")] <- NA_character_
meta[[gcol]] <- tmp

# Align meta to matrix
meta <- meta[match(colnames(M), meta$SampleID), , drop = FALSE]
if (any(is.na(meta$SampleID))) stop("Metadata template is missing some SampleIDs present in matrix.")
g_factor <- factor(meta[[gcol]], levels = c("HER2neg","HER2pos"))
if (any(is.na(g_factor))) stop("Metadata template has missing/invalid groups.")
if (min(table(g_factor)) < min_per_group) stop("HER2 groups too small after meta template: ", paste(table(g_factor), collapse=" / "))

# ---------------- 5) CLEAN/IMPUTE/FILTER ----------------
M <- robust_impute_and_filter(M)
stopifnot(nrow(M) >= 2)

# ---------------- 6) LIMMA DEP ----------------
design <- model.matrix(~ 0 + g_factor); colnames(design) <- c("HER2neg","HER2pos")
fit  <- lmFit(M, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(HER2pos - HER2neg, levels = design)),
               trend = TRUE, robust = TRUE)

tt <- topTable(fit2, number = Inf, sort.by = "P")
tt$protein <- rownames(tt)
tt$log2FC_HER2pos_vs_HER2neg <- tt$logFC
attr(tt, "Amean") <- fit$Amean

# ---------------- 7) WRITER (Krug-style artifacts) ----------------
write_outputs <- function(res, tag, E_used, g_used){
  out_dir <- out_dir_tag(tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # meta USED + contrast USED
  meta_used_tsv     <- file.path(out_dir, paste0(out_prefix, "_", tag, "_meta_USED.tsv"))
  contrast_used_txt <- file.path(out_dir, paste0(out_prefix, "_", tag, "_contrast_USED.txt"))
  
  meta_used <- data.table(
    SampleID = colnames(E_used),
    Group    = as.character(g_used)
  )
  # optional extra context if present in template
  if ("HER2_raw" %in% names(meta)) meta_used[, HER2_raw := meta$HER2_raw]
  if ("QC_raw"   %in% names(meta)) meta_used[, QC_raw   := meta$QC_raw]
  meta_used[, QC_filter := ifelse(USE_QC_PASS, "pass_only", "not_applied")]
  
  fwrite(meta_used, meta_used_tsv, sep = "\t")
  
  writeLines(c(
    "Mode: HER2",
    "Contrast: HER2pos - HER2neg",
    paste0("Definition tag: ", tag),
    paste0("HER2 definition: Supp01 HER2 Status (clinical)"),
    paste0("QC filter: ", ifelse(USE_QC_PASS, "pass_only", "not_applied")),
    paste0("n_pos: ", sum(g_used == "HER2pos")),
    paste0("n_neg: ", sum(g_used == "HER2neg"))
  ), contrast_used_txt)
  
  # decorate results
  res$protein         <- toupper(res$protein)
  res$cohort_id       <- paste0(out_prefix, "_", tag)
  res$HER2_definition <- "Supp01 HER2 Status (clinical)"
  res$QC_filter       <- ifelse(USE_QC_PASS, "pass_only", "not_applied")
  res$feature_id_source <- mx$feature_col
  res$n_pos <- as.integer(sum(g_used == "HER2pos"))
  res$n_neg <- as.integer(sum(g_used == "HER2neg"))
  
  # full
  full_path <- file.path(out_dir, "DE_full.tsv")
  fwrite(as.data.table(res), full_path, sep = "\t")
  
  # filtered tables (Krug-style names)
  sig <- res[is.finite(res$adj.P.Val) & res$adj.P.Val < FDR_MAX,
             c("protein","log2FC_HER2pos_vs_HER2neg","AveExpr","t","P.Value","adj.P.Val","B",
               "cohort_id","HER2_definition","QC_filter","feature_id_source","n_pos","n_neg")]
  
  f_fdr <- file.path(out_dir, paste0("DEP_HER2pos_vs_HER2neg_", tag, "_FDRlt0.05.tsv"))
  fwrite(sig, f_fdr, sep = "\t")
  
  sig_fc2 <- sig[is.finite(sig$log2FC_HER2pos_vs_HER2neg) &
                   abs(sig$log2FC_HER2pos_vs_HER2neg) >= FC_MIN_LOG2, ]
  f_fc2 <- file.path(out_dir, paste0("DEP_HER2pos_vs_HER2neg_", tag, "_FDRlt0.05_absLog2FCge1.tsv"))
  fwrite(sig_fc2, f_fc2, sep = "\t")
  
  # canonical copies (Krug-style)
  canon_fdr <- file.path(out_dir, paste0(out_prefix, "_", tag, "_DE_filtered_FDRlt0.05.tsv"))
  canon_fc2 <- file.path(out_dir, paste0(out_prefix, "_", tag, "_DE_filtered_FDRlt0.05_absLog2FCge1.tsv"))
  file.copy(f_fdr, canon_fdr, overwrite = TRUE)
  file.copy(f_fc2, canon_fc2, overwrite = TRUE)
  
  # volcano (STANDARDIZED: -log10(FDR))
  safe_fdr <- res$adj.P.Val
  safe_fdr[!is.finite(safe_fdr)] <- NA_real_
  safe_fdr[safe_fdr == 0] <- .Machine$double.xmin
  
  volc_path <- file.path(out_dir, "volcano.pdf")
  pdf(volc_path, 7.5, 6)
  par(mar = c(5,5,4.5,2))
  
  x <- res$log2FC_HER2pos_vs_HER2neg
  y <- -log10(safe_fdr)
  
  is_sig <- is.finite(res$adj.P.Val) & res$adj.P.Val < FDR_MAX
  is_fc  <- is.finite(x) & abs(x) >= FC_MIN_LOG2
  
  plot(x, y,
       xlab = "log2FC (HER2pos - HER2neg)",
       ylab = "-log10(FDR)",
       main = paste0("Volcano — ", tag),
       pch = 21,
       bg  = ifelse(is_sig, "black", "white"),
       col = "black")
  abline(h = -log10(FDR_MAX), v = c(-FC_MIN_LOG2, FC_MIN_LOG2), lty = 2, col = "grey50")
  grid(col = "grey90")
  
  labs <- res$protein
  is_bad <- grepl("^UNMAPPED_", labs) | is.na(labs) | labs == ""
  labs[is_bad] <- NA
  
  force <- c("ERBB2","GRB7","STARD3","MIEN1","ERBB4")
  idx_force <- which(!is.na(labs) & labs %in% force)
  
  idx_core <- which(is_sig & is_fc & !is_bad)
  if (length(idx_core) > 5) idx_core <- idx_core[order(res$adj.P.Val[idx_core])][1:5]
  
  idx_top <- order(res$adj.P.Val)
  idx_top <- idx_top[!idx_top %in% c(idx_force, idx_core)]
  idx_top <- idx_top[!is_bad[idx_top]]
  idx_top <- head(idx_top, 12)
  
  idx_lab <- unique(c(idx_force, idx_core, idx_top))
  idx_lab <- idx_lab[!is.na(labs[idx_lab])]
  
  if (length(idx_lab)) {
    dx <- rep(c(-0.08, 0.08), length.out = length(idx_lab))
    dy <- rep(c(0.12, 0.18), length.out = length(idx_lab))
    text(x[idx_lab] + dx, y[idx_lab] + dy, labels = labs[idx_lab],
         cex = 0.8, pos = ifelse(dx > 0, 4, 2), xpd = NA)
    points(x[idx_lab], y[idx_lab], pch = 21,
           bg = ifelse(is_sig[idx_lab] & is_fc[idx_lab], "black", "white"),
           col = "black")
  }
  dev.off()
  
  # MA
  ma_path <- file.path(out_dir, "MA.pdf")
  pdf(ma_path, 7.5, 6)
  plot(attr(res, "Amean"), res$log2FC_HER2pos_vs_HER2neg,
       xlab = "AveExpr", ylab = "log2FC",
       main = paste0("MA — ", tag),
       pch = 21, bg = "white", col = "black")
  abline(h = 0, lty = 2, col = "grey50")
  grid(col = "grey90")
  dev.off()
  
  # QC
  qc_box(E_used, g_used, tag, "ERBB2")
  qc_box(E_used, g_used, tag, "GRB7")
  
  # RUN_LOG
  erbb2_line <- {
    i <- which(res$protein == "ERBB2")[1]
    if (!is.na(i)) sprintf("ERBB2: log2FC=%.3f  FDR=%.3g",
                           res$log2FC_HER2pos_vs_HER2neg[i], res$adj.P.Val[i])
    else "ERBB2: NA"
  }
  grb7_line <- {
    i <- which(res$protein == "GRB7")[1]
    if (!is.na(i)) sprintf("GRB7:  log2FC=%.3f  FDR=%.3g",
                           res$log2FC_HER2pos_vs_HER2neg[i], res$adj.P.Val[i])
    else "GRB7:  NA"
  }
  
  log_lines <- c(
    sprintf("RUN LOG — %s_%s", out_prefix, tag),
    sprintf("Supp01: %s", basename(supp1_xlsx)),
    sprintf("Supp03: %s", basename(supp3_xlsx)),
    sprintf("Feature column used: %s", mx$feature_col),
    sprintf("QC filter: %s", ifelse(USE_QC_PASS, "pass_only", "not_applied")),
    sprintf("Samples kept: %d (HER2pos=%d, HER2neg=%d)", length(g_used), sum(g_used=="HER2pos"), sum(g_used=="HER2neg")),
    sprintf("FDR<0.05: %d", nrow(sig)),
    sprintf("FDR<0.05 & |log2FC|>=1: %d", nrow(sig_fc2)),
    "Meta/contrast:",
    paste0("  - ", meta_template_csv),
    paste0("  - ", meta_used_tsv),
    paste0("  - ", contrast_used_txt),
    "Outputs:",
    paste0("  - ", full_path),
    paste0("  - ", f_fdr),
    paste0("  - ", f_fc2),
    paste0("  - ", canon_fdr),
    paste0("  - ", canon_fc2),
    paste0("  - ", volc_path),
    paste0("  - ", ma_path),
    paste0("  - ", file.path(out_dir, "QC_ERBB2_boxplot.pdf")),
    paste0("  - ", file.path(out_dir, "QC_GRB7_boxplot.pdf")),
    erbb2_line,
    grb7_line
  )
  
  writeLines(log_lines, file.path(out_dir, "RUN_LOG.txt"))
  invisible(list(n_fdr = nrow(sig), n_fc2 = nrow(sig_fc2)))
}

# ---------------- RUN ----------------
stats <- write_outputs(tt, TAG, M, as.character(g_factor))

cat(sprintf("[DEP %s] FDR<0.05=%d | FDR+FC2=%d\n", TAG, stats$n_fdr, stats$n_fc2))
# ============================ END ============================
