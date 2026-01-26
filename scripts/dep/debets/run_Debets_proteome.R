# ================== run_Debets_proteome.R — FINAL (FORCE pCR; Debets HER2+ cohort) ==================
# Debets proteome is NOT HER2pos vs HER2neg — it’s response (pCR vs No_pCR) within HER2+ tumors.
# Produces:
#   - meta_USED.tsv + contrast_USED.txt (what limma actually used)
#   - DE_full.tsv + filtered (FDR + FDR+FC2) + canonical copies
#   - volcano + ERBB2/GRB7 QC + RUN_LOG

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(matrixStats)
  library(limma)
  library(stringr)
})

# ================== CONFIG OVERRIDES ==================
FORCE_PCR_MODE <- TRUE
# ======================================================

# --- locate this script's folder ---
here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) && file.exists(f)) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(".")
})
base <- here
message("Base: ", base)

# --- file paths ---
PROT_X  <- file.path(base, "table_s2_Protein_Quant.xlsx")
META_C  <- file.path(base, "Debets2023_protein_meta_template.csv")  # <-- uses the filled one
OUT_TAG <- "Debets2023_protein"
OUT_DIR <- file.path(base, "aim2_outputs", OUT_TAG)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- params ---
NONMISS_FRAC   <- 0.60  # fraction of non-missing per protein across all samples
MIN_PER_GROUP  <- 2     # require >= this many non-NA per group (after filtering)
MERGE_nPCR_INTO_No <- TRUE
FC_MIN_LOG2    <- 1
FDR_MAX        <- 0.05

stopifnot(file.exists(PROT_X), file.exists(META_C))

# ---------------- helpers ----------------
pool_token <- function(s) {
  s <- gsub("[[:cntrl:]]", "", s)
  s <- trimws(as.character(s))
  m <- regexpr("Pool\\d+_[0-9]{3}[NC]?", s, perl = TRUE, ignore.case = TRUE)
  out <- ifelse(m > 0, substr(s, m, m + attr(m, "match.length") - 1), NA_character_)
  out <- sub("^pool", "Pool", out, ignore.case = TRUE)
  out <- sub("([0-9]{3})([nc])$", "\\1\\U\\2", out, perl = TRUE)
  out
}

read_quant_xlsx <- function(path) {
  peek <- suppressMessages(readxl::read_excel(path, col_names = FALSE, n_max = 10))
  peek <- as.data.frame(peek)
  pool_counts <- sapply(seq_len(nrow(peek)), function(r)
    sum(!is.na(pool_token(unlist(peek[r, ], use.names = FALSE)))))
  hdr_row <- which.max(pool_counts)
  if (pool_counts[hdr_row] < 3) stop("Can't find a header row with Pool labels in: ", basename(path))
  
  full <- suppressMessages(readxl::read_excel(path, col_names = FALSE))
  full <- as.data.frame(full)
  
  colnames(full) <- make.names(unlist(full[hdr_row, ], use.names = FALSE), unique = TRUE)
  full <- full[(hdr_row + 1):nrow(full), , drop = FALSE]
  full <- full[, colSums(!is.na(full)) > 0, drop = FALSE]
  full
}

to_numeric_df <- function(df) {
  for (j in seq_len(ncol(df))) {
    x <- as.character(df[[j]])
    x[x %in% c("", "NA", "NaN", "na", "N/A")] <- NA_character_
    suppressWarnings(df[[j]] <- as.numeric(x))
  }
  df
}

norm_str <- function(x) {
  x <- trimws(as.character(x)); x <- tolower(x)
  x <- gsub("\\s+", "", x); x <- gsub("[^a-z0-9+_]", "", x)
  x
}

map_outcome <- function(v) {
  raw <- norm_str(v)
  lab <- rep(NA_character_, length(raw))
  
  lab[raw %in% c("pcr","pathcr","completepcr","completepathcr","cpr")] <- "pCR"
  lab[raw %in% c("nopcr","nonpcr","noresponse","nonresponse","nonresponder","residualdisease","rd","no_pcr")] <- "No_pCR"
  
  lab[grepl("npcr|near", raw)] <- if (MERGE_nPCR_INTO_No) "No_pCR" else "npCR"
  
  lab[is.na(lab) & grepl("pcr", raw) & !grepl("(^|[^a-z])(no|non|not)", raw)] <- "pCR"
  lab[is.na(lab) & grepl("pcr", raw) &  grepl("(no|non|not|residual|rd)", raw)] <- "No_pCR"
  
  lab
}

safe_sd <- function(x) {
  s <- suppressWarnings(sd(x, na.rm = TRUE))
  if (!is.finite(s)) 0 else s
}

# ---------------- 1) load proteome table ----------------
raw_df <- read_quant_xlsx(PROT_X)

nm <- names(raw_df)
if ("Gene.Name" %in% nm) {
  id_colname <- "Gene.Name"
} else if ("Protein.ID" %in% nm) {
  id_colname <- "Protein.ID"
} else {
  is_texty <- sapply(raw_df, function(z) {
    zc <- suppressWarnings(as.numeric(as.character(z)))
    mean(is.na(zc)) > 0.5
  })
  id_col <- which(is_texty)[1]
  if (!length(id_col)) stop("No obvious ID column.")
  id_colname <- names(raw_df)[id_col]
}

row_ids <- make.unique(as.character(raw_df[[id_colname]]))
raw_df[[id_colname]] <- NULL

tok <- pool_token(names(raw_df))
meas <- raw_df[, !is.na(tok), drop = FALSE]
names(meas) <- tok[!is.na(tok)]
meas <- to_numeric_df(meas)
meas <- meas[, colSums(!is.na(meas)) > 0, drop = FALSE]
if (ncol(meas) == 0) stop("No usable measurement columns.")

X <- as.matrix(meas)
storage.mode(X) <- "double"
rownames(X) <- row_ids

# ---------------- 2) load metadata ----------------
meta <- fread(META_C)
stopifnot("SampleID" %in% names(meta), "Outcome" %in% names(meta))

meta$SampleID <- pool_token(meta$SampleID)

if (!("Plex" %in% names(meta)) || all(is.na(meta$Plex) | !nzchar(meta$Plex))) {
  meta$Plex <- sub("^(Pool\\d+).*", "\\1", meta$SampleID)
}

# ---------------- 3) choose analysis mode (FORCED pCR) ----------------
grp_OUT <- map_outcome(meta$Outcome)

if (isTRUE(FORCE_PCR_MODE)) {
  ANALYSIS <- "pCR"
  meta$Group <- ifelse(grp_OUT %in% c("pCR","No_pCR"), grp_OUT, NA_character_)
  CONTRAST_LABEL <- "pCR_vs_No_pCR"
  message("FORCING pCR mode for Debets dataset (within HER2+ tumors).")
} else {
  stop("This script is intended for Debets and should be forced into pCR mode.")
}

if (length(unique(na.omit(meta$Group))) < 2) {
  stop("pCR mode selected but <2 outcome groups present after mapping. Check Outcome column values.")
}

# ---------------- 4) align samples ----------------
common <- intersect(colnames(X), meta$SampleID)
if (length(common) < 6) stop("Not enough overlap between data and metadata.")

X <- X[, common, drop = FALSE]
meta <- meta[match(common, meta$SampleID), , drop = FALSE]

keep <- meta$Group %in% c("pCR","No_pCR")
X <- X[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]

if (ncol(X) < 6) stop("Need >=6 samples after filtering.")

message("Group counts:"); print(table(meta$Group, useNA = "no"))

# ================== EXPORT LIMMA INPUTS (Nick-ready) ==================
meta_export <- copy(meta)
keep_cols <- intersect(
  c("SampleID","Group","Outcome","HER2_IHC","HER2_RATIO","Plex","Bridge",
    "Patient","Histology","Tumour_grade","ER_pct","PR_pct","Tumour_pct"),
  names(meta_export)
)
meta_export <- meta_export[, ..keep_cols]

meta_outfile <- file.path(OUT_DIR, paste0(OUT_TAG, "_meta_USED.tsv"))
fwrite(meta_export, meta_outfile, sep = "\t")

contrast_used <- file.path(OUT_DIR, paste0(OUT_TAG, "_contrast_USED.txt"))
contrast_string <- "GrouppCR - GroupNo_pCR"
writeLines(c(
  paste0("Mode: ", ANALYSIS),
  paste0("Contrast: ", contrast_string),
  paste0("Contrast label: ", CONTRAST_LABEL),
  paste0("n_pCR: ", sum(meta$Group == "pCR")),
  paste0("n_No_pCR: ", sum(meta$Group == "No_pCR"))
), contrast_used)

# ---------------- 5) within-plex centering + batch removal ----------------
X1 <- X
for (p in unique(meta$Plex)) {
  idx <- which(meta$Plex == p)
  cm  <- matrixStats::colMedians(X1[, idx, drop = FALSE], na.rm = TRUE)
  X1[, idx] <- sweep(X1[, idx, drop = FALSE], 2, cm, "-")
}

if (length(unique(meta$Plex)) > 1) {
  des_batch <- model.matrix(~ Group, data = meta)
  X2 <- limma::removeBatchEffect(X1, batch = factor(meta$Plex), design = des_batch)
} else {
  X2 <- X1
}

# ---------------- 6) filter proteins ----------------
meta$Group <- factor(meta$Group, levels = c("No_pCR", "pCR"))
grp <- meta$Group

lv <- levels(grp)
ok_each_group <-
  rowSums(!is.na(X2[, grp == lv[1], drop = FALSE])) >= MIN_PER_GROUP &
  rowSums(!is.na(X2[, grp == lv[2], drop = FALSE])) >= MIN_PER_GROUP

keep_rows <- rowMeans(!is.na(X2)) >= NONMISS_FRAC & ok_each_group
Xf <- X2[keep_rows, , drop = FALSE]
Xf <- Xf[apply(Xf, 1, safe_sd) > 0, , drop = FALSE]
if (nrow(Xf) == 0) stop("No proteins left after filtering. Relax NONMISS_FRAC or MIN_PER_GROUP.")

# ---------------- 7) limma DE ----------------
tab <- table(meta$Plex, meta$Group)
plex_counts <- rowSums(tab)
has_both <- rowSums(tab > 0) == 2
use_plex <- (nrow(tab) > 1) && all(plex_counts >= 2) && any(has_both)

if (use_plex) {
  design <- model.matrix(~ 0 + Group + Plex, data = meta)
} else {
  design <- model.matrix(~ 0 + Group, data = meta)
}
colnames(design) <- make.names(colnames(design))

fit  <- lmFit(Xf, design)
ct   <- makeContrasts(Contrast = GrouppCR - GroupNo_pCR, levels = design)
fit2 <- eBayes(contrasts.fit(fit, ct), trend = TRUE, robust = TRUE)

res <- topTable(fit2, number = Inf, sort.by = "P")
res$feature_id <- rownames(res)
res$log2FC_pCR_vs_No_pCR <- res$logFC

# ---------------- 8) outputs (match Krug/Rajkumar style) ----------------
full_cols <- c("feature_id","log2FC_pCR_vs_No_pCR","AveExpr","t","P.Value","adj.P.Val","B")
res_full <- res[, c("feature_id","log2FC_pCR_vs_No_pCR","AveExpr","t","P.Value","adj.P.Val","B"), drop = FALSE]
fwrite(res_full, file.path(OUT_DIR, "DE_full.tsv"), sep = "\t")

sig <- res[is.finite(res$adj.P.Val) & res$adj.P.Val < FDR_MAX, ]
sig_out <- data.table(
  feature_id = sig$feature_id,
  log2FC_pCR_vs_No_pCR = sig$log2FC_pCR_vs_No_pCR,
  AveExpr = sig$AveExpr, t = sig$t,
  P.Value = sig$P.Value, adj.P.Val = sig$adj.P.Val, B = sig$B,
  cohort_id = OUT_TAG,
  outcome_definition = "Debets clinical outcome",
  n_pCR = sum(meta$Group == "pCR"),
  n_No_pCR = sum(meta$Group == "No_pCR")
)

f_fdr <- file.path(OUT_DIR, paste0("DEP_", CONTRAST_LABEL, "_FDRlt0.05.tsv"))
fwrite(sig_out, f_fdr, sep = "\t")

sig_fc2 <- sig_out[is.finite(log2FC_pCR_vs_No_pCR) & abs(log2FC_pCR_vs_No_pCR) >= FC_MIN_LOG2]
f_fc2 <- file.path(OUT_DIR, paste0("DEP_", CONTRAST_LABEL, "_FDRlt0.05_absLog2FCge1.tsv"))
fwrite(sig_fc2, f_fc2, sep = "\t")

# canonical copies
canon_fdr <- file.path(OUT_DIR, paste0(OUT_TAG, "_DE_filtered_FDRlt0.05.tsv"))
canon_fc2 <- file.path(OUT_DIR, paste0(OUT_TAG, "_DE_filtered_FDRlt0.05_absLog2FCge1.tsv"))
file.copy(f_fdr, canon_fdr, overwrite = TRUE)
file.copy(f_fc2, canon_fc2, overwrite = TRUE)

# ---------------- 9) QC ERBB2 / GRB7 ----------------
qc_box <- function(E, g_factor, gene) {
  gene <- toupper(gene)
  if (!(gene %in% rownames(E))) return(invisible(FALSE))
  
  out_pdf <- file.path(OUT_DIR, paste0("QC_", gene, "_boxplot.pdf"))
  vals <- as.numeric(E[gene, ])
  v_pos <- vals[g_factor == "pCR"]
  v_neg <- vals[g_factor == "No_pCR"]
  m_pos <- mean(v_pos, na.rm = TRUE); s_pos <- sd(v_pos, na.rm = TRUE)
  m_neg <- mean(v_neg, na.rm = TRUE); s_neg <- sd(v_neg, na.rm = TRUE)
  log2fc <- m_pos - m_neg
  ttst <- try(stats::t.test(v_pos, v_neg, var.equal = FALSE), silent = TRUE)
  pval <- if (inherits(ttst, "try-error")) NA_real_ else ttst$p.value
  
  pdf(out_pdf, 6.8, 5.6)
  par(mar = c(5,5,4.5,2))
  boxplot(split(vals, g_factor),
          main = sprintf("%s QC (Debets proteome)", gene),
          ylab = "Expression (centered / batch-adjusted)", xlab = "Outcome group",
          col = c("grey85","grey85"), border = "grey25")
  stripchart(split(vals, g_factor), vertical = TRUE, method = "jitter",
             pch = 21, bg = "white", col = "black", add = TRUE)
  grid(col = "grey90")
  mtext(sprintf("n(pCR)=%d  n(No_pCR)=%d   mean±SD(pCR)=%.2f±%.2f   mean±SD(No_pCR)=%.2f±%.2f   log2FC=%.2f   t-test P=%.3g",
                length(v_pos), length(v_neg), m_pos, s_pos, m_neg, s_neg, log2fc, pval),
        side = 3, line = 0.6, cex = 0.85)
  dev.off()
  TRUE
}

qc_box(Xf, meta$Group, "ERBB2")
qc_box(Xf, meta$Group, "GRB7")

# ---------------- 10) Volcano ----------------
p05 <- -log10(FDR_MAX)

pdf(file.path(OUT_DIR, "volcano.pdf"), 8.5, 6.5)
par(mar = c(5,5,4.5,2))

x <- res$logFC
y <- -log10(res$P.Value)
labs <- res$feature_id

is_sig_fdr <- res$adj.P.Val < FDR_MAX
is_fc1     <- abs(res$logFC) >= FC_MIN_LOG2

plot(x, y,
     pch = 21,
     bg  = ifelse(is_sig_fdr, "black", "white"),
     col = "black",
     xlab = "log2FC (pCR - No_pCR)",
     ylab = "-log10(P)",
     main = OUT_TAG)

abline(v = c(-FC_MIN_LOG2, FC_MIN_LOG2), lty = 2, col = "grey55")
abline(h = p05, lty = 2, col = "grey55")
grid(col = "grey90")

is_bad_label <- is.na(labs) | labs == ""
labs[is_bad_label] <- NA

force_genes <- c("ERBB2", "GRB7", "STARD3", "MIEN1")
idx_force <- which(!is.na(labs) & labs %in% force_genes)

idx_core <- which(is_sig_fdr & is_fc1 & !is_bad_label)
max_core_labels <- 5
if (length(idx_core) > max_core_labels) idx_core <- idx_core[order(res$P.Value[idx_core])][1:max_core_labels]

topN_extra <- 12
idx_topP <- order(res$P.Value)
idx_topP <- idx_topP[!idx_topP %in% c(idx_force, idx_core)]
idx_topP <- idx_topP[!is_bad_label[idx_topP]]
idx_topP <- head(idx_topP, topN_extra)

idx_label <- unique(c(idx_force, idx_core, idx_topP))
idx_label <- idx_label[!is.na(labs[idx_label])]

dx <- rep(c(-0.08, 0.08), length.out = length(idx_label))
dy <- rep(c(0.12, 0.18), length.out = length(idx_label))

text(x[idx_label] + dx, y[idx_label] + dy,
     labels = labs[idx_label],
     cex = 0.75,
     pos = ifelse(dx > 0, 4, 2),
     xpd = NA)

points(x[idx_label], y[idx_label],
       pch = 21,
       bg = ifelse(is_sig_fdr[idx_label] & is_fc1[idx_label], "black", "white"),
       col = "black")

dev.off()

# ---------------- 11) RUN LOG ----------------
erbb2_line <- {
  i <- which(res$feature_id == "ERBB2")[1]
  if (!is.na(i)) sprintf("ERBB2: log2FC=%.3f  FDR=%.3g", res$logFC[i], res$adj.P.Val[i]) else "ERBB2: NA"
}
grb7_line <- {
  i <- which(res$feature_id == "GRB7")[1]
  if (!is.na(i)) sprintf("GRB7:  log2FC=%.3f  FDR=%.3g", res$logFC[i], res$adj.P.Val[i]) else "GRB7: NA"
}

log_lines <- c(
  sprintf("RUN LOG — %s", OUT_TAG),
  sprintf("Proteome file: %s", PROT_X),
  sprintf("Meta template used: %s", META_C),
  sprintf("Mode: %s (contrast: %s)", ANALYSIS, CONTRAST_LABEL),
  sprintf("Samples kept: %d (pCR=%d, No_pCR=%d)", ncol(X), sum(meta$Group=="pCR"), sum(meta$Group=="No_pCR")),
  sprintf("Plex adjustment used: %s", ifelse(use_plex, "Group + Plex", "Group only")),
  sprintf("Filtering: NONMISS_FRAC=%.2f  MIN_PER_GROUP=%d", NONMISS_FRAC, MIN_PER_GROUP),
  sprintf("FDR<%.2f: %d", FDR_MAX, nrow(sig_out)),
  sprintf("FDR<%.2f & |log2FC|>=%.1f: %d", FDR_MAX, FC_MIN_LOG2, nrow(sig_fc2)),
  "Meta artifacts:",
  paste0("  - ", meta_outfile),
  paste0("  - ", contrast_used),
  "Outputs:",
  paste0("  - ", file.path(OUT_DIR, "DE_full.tsv")),
  paste0("  - ", f_fdr),
  paste0("  - ", f_fc2),
  paste0("  - ", canon_fdr),
  paste0("  - ", canon_fc2),
  paste0("  - ", file.path(OUT_DIR, "volcano.pdf")),
  erbb2_line,
  grb7_line
)
writeLines(log_lines, file.path(OUT_DIR, "RUN_LOG.txt"))

message("Done. Wrote FDR (n=", nrow(sig_out), ") and FC2 (n=", nrow(sig_fc2), ").")
message("Complete → ", OUT_DIR)
# ================== end ==================
