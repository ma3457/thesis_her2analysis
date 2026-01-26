# =========== run_GSE81538_RNA.R (Brueffer) ===============

suppressPackageStartupMessages({
  library(data.table); library(limma); library(pheatmap); library(stringr)
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
expr_file <- file.path(script_dir, "GSE81538_gene_expression_405_transformed.csv")
lab_file  <- file.path(script_dir, "GSE81538_prep", "HER2_mapping_matched.csv")
out_tag   <- "GSE81538"
out_dir   <- file.path(script_dir, "aim2_outputs", out_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- meta USED + contrast USED paths ---
meta_used_tsv     <- file.path(out_dir,  "GSE81538_RNA_meta_USED.tsv")
contrast_used_txt <- file.path(out_dir,  "GSE81538_RNA_contrast_USED.txt")

wrap_title <- function(s, width = 60) paste(strwrap(s, width = width), collapse = "\n")
plot_with_margins <- function(expr) { op <- par(mar = c(5,5,5,2)); on.exit(par(op), add = TRUE); force(expr) }

norm_ids <- function(v) {
  v <- as.character(v)
  v <- sub("^X", "", v)
  v <- sub("\\.GPL\\d+$", "", v)
  v <- gsub("[^A-Za-z0-9]", "", v)
  toupper(v)
}

# ---------------- read expression ----------------
if (!file.exists(expr_file)) stop("Expression matrix not found: ", expr_file)
x <- fread(expr_file, check.names = FALSE)
if (is.na(names(x)[1]) || names(x)[1] == "") setnames(x, 1, "gene_symbol")
if (!names(x)[1] %in% c("gene_symbol","Gene","feature_id","ID","Feature")) setnames(x, 1, "gene_symbol")
x <- x[!duplicated(x$gene_symbol), ]

# coerce numeric sample columns safely
samp_cols <- setdiff(names(x), "gene_symbol")
for (nm in samp_cols) {
  if (!is.numeric(x[[nm]])) {
    v <- as.character(x[[nm]]); v[is.na(v)] <- ""
    suppressWarnings(num <- as.numeric(v))
    if (sum(nchar(v) > 0 & is.na(num)) == 0) x[[nm]] <- num
  }
}
E <- as.matrix(x[, ..samp_cols]); rownames(E) <- x$gene_symbol
storage.mode(E) <- "numeric"

# ---------------- read HER2 labels (Brueffer) ----------------
if (!file.exists(lab_file)) stop("Label file not found: ", lab_file)
lab <- fread(lab_file, check.names = FALSE)

# Pick an ID column and rename to sample_id
id_col <- intersect(names(lab), c("sample_id","SampleID","SampleID_raw","Sample","sample","geo_accession","ID"))
if (!length(id_col)) id_col <- names(lab)[1]
setnames(lab, id_col[1], "sample_id")

# HER2 status -> group (Positive/Negative), robust to different column names
her2_col <- intersect(names(lab), c("HER2_status","her2_status","HER2","her2","group"))
if (!length(her2_col)) {
  stop("Could not find a HER2 status column in HER2_mapping_matched.csv.\nColumns: ",
       paste(names(lab), collapse = ", "))
}
if (her2_col[1] != "group") setnames(lab, her2_col[1], "group")

lab$sample_id <- str_trim(as.character(lab$sample_id))
lab$group     <- str_trim(as.character(lab$group))

# Map group to Positive/Negative if not already
lab$group <- ifelse(
  grepl("pos|3\\+|amp|amplif|overexpr|her2[-_ ]?high", lab$group, ignore.case = TRUE),
  "Positive",
  ifelse(
    grepl("neg|0|1\\+|her2[-_ ]?low", lab$group, ignore.case = TRUE),
    "Negative",
    NA_character_
  )
)

stopifnot("sample_id" %in% names(lab))
stopifnot("group" %in% names(lab))

# ---------------- align by normalized IDs ----------------
expr_ids_norm <- norm_ids(colnames(E))
lab_ids_norm  <- norm_ids(lab$sample_id)

map_dt  <- data.table(sample_id = lab$sample_id, sample_norm = lab_ids_norm)
expr_dt <- data.table(expr_col  = colnames(E),   sample_norm = expr_ids_norm)

m <- merge(map_dt, expr_dt, by = "sample_norm")
if (nrow(m)) {
  m <- m[, .SD[.N == 1], by = sample_norm]
  lab <- merge(lab, m[, .(sample_id, expr_col)], by = "sample_id", all.x = TRUE)
  lab$sample_id <- ifelse(!is.na(lab$expr_col), lab$expr_col, lab$sample_id)
  lab$expr_col <- NULL
}

keep <- intersect(colnames(E), lab$sample_id)
if (!length(keep)) stop("No overlap between expression columns and labels$sample_id")

# finalize alignment + drop NA labels / non-Pos/Neg
E   <- E[, keep, drop = FALSE]
lab <- lab[match(colnames(E), lab$sample_id), , drop = FALSE]
ok  <- !is.na(lab$group) & lab$group %in% c("Positive","Negative")
E   <- E[, ok, drop = FALSE]
lab <- lab[ok, , drop = FALSE]

# groups
grp <- factor(ifelse(lab$group == "Positive", "HER2pos","HER2neg"),
              levels = c("HER2neg","HER2pos"))
n_pos <- sum(grp == "HER2pos"); n_neg <- sum(grp == "HER2neg")
if (n_pos < 3 || n_neg < 3) stop("Groups too small (need ≥3 per arm).")

# ---------------- write META USED + CONTRAST USED ----------------
meta_used <- data.table(
  SampleID  = colnames(E),
  Group     = as.character(grp),
  Group_raw = lab$group
)
fwrite(meta_used, meta_used_tsv, sep = "\t", quote = FALSE, na = "NA")
writeLines(c(
  "Mode: HER2",
  "Contrast: HER2pos - HER2neg",
  paste0("Cohort: ", out_tag),
  paste0("n_pos: ", n_pos),
  paste0("n_neg: ", n_neg)
), contrast_used_txt)
message("Wrote meta USED: ", meta_used_tsv)
message("Wrote contrast USED: ", contrast_used_txt)

# quick log2 if clearly raw
vals <- as.numeric(E[is.finite(E)])
if (quantile(vals, 0.99, na.rm = TRUE) > 30 || max(vals, na.rm = TRUE) > 1000) {
  E <- log2(E + 1)
}

# drop zero-variance
zv <- apply(E, 1, sd, na.rm = TRUE) == 0
if (any(zv)) E <- E[!zv, , drop = FALSE]

# ---------------- LIMMA (HER2pos - HER2neg) ----------------
design <- model.matrix(~ 0 + grp); colnames(design) <- c("HER2neg","HER2pos")
fit  <- lmFit(E, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(HER2pos - HER2neg, levels = design)),
               trend = TRUE, robust = TRUE)
tt <- topTable(fit2, number = Inf, sort.by = "P")
tt$feature_id <- rownames(tt)

# ---------------- standard outputs (Wolf-style) ----------------
all_path <- file.path(out_dir, "ALL_HER2pos_vs_HER2neg_DE.tsv")
fwrite(tt, all_path, sep = "\t")

tt_fdr <- tt[tt$adj.P.Val < 0.05, ]
canon_fdr <- file.path(out_dir, "GSE81538_DE_filtered_FDRlt0.05.tsv")
fwrite(tt_fdr, canon_fdr, sep = "\t")

tabf <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1.0, ]
canon_fc2 <- file.path(out_dir, "GSE81538_DE_filtered_FDRlt0.05_absLog2FCge1.tsv")
fwrite(tabf, canon_fc2, sep = "\t")

filt_path <- file.path(out_dir, "FILT_HER2pos_vs_HER2neg_DE_FDR0.05_LFC1.tsv")
fwrite(tabf[, c("feature_id","logFC","adj.P.Val")], filt_path, sep = "\t")

if (nrow(tabf) >= 3) {
  keep_genes <- intersect(tabf$feature_id, rownames(E))
  if (length(keep_genes) >= 3) {
    M <- E[keep_genes, , drop = FALSE]
    cm <- cor(t(M), use = "pairwise.complete.obs")
    png(file.path(out_dir, "CORR_HER2pos_vs_HER2neg_filtered.png"),
        width = 1400, height = 1200, res = 180)
    pheatmap(cm, main = "Correlation — HER2pos vs HER2neg (FDR<0.05 & |log2FC|≥1)")
    dev.off()
  }
}

# ---------------- quick QC plots (ERBB2 / GRB7) ----------------
qc_genes <- c("ERBB2","GRB7")
for (g in qc_genes) {
  if (g %in% rownames(E)) {
    ttl <- wrap_title(sprintf("%s %s | Brueffer (n+=%d, n-=%d)", out_tag, g, n_pos, n_neg), 55)
    fn <- file.path(out_dir, sprintf("qc_%s_boxplot.png", g))
    png(fn, width = 1400, height = 1100, res = 150)
    plot_with_margins({
      boxplot(split(as.numeric(E[g, ]), grp),
              main = ttl, ylab = "Expression", xlab = "HER2 group",
              col="grey90", border="grey30", cex.main = 0.95)
    })
    dev.off()
  }
}

# ---------------- Volcano (FDR) + MA (PNG) ----------------
sig_flag <- tt$adj.P.Val < 0.05
p05_fdr  <- -log10(0.05)

plot_volcano <- function() {
  ttl <- wrap_title(sprintf("%s Volcano (Brueffer, n+=%d, n-=%d)", out_tag, n_pos, n_neg), 55)
  
  y <- -log10(pmax(tt$adj.P.Val, 1e-300))  # FDR
  plot(tt$logFC, y,
       xlab="log2FC (pos vs neg)",
       ylab=expression(-log[10]("FDR")),
       main=ttl,
       pch=21,
       bg=ifelse(sig_flag, "black", "white"),
       col="black",
       cex.main=0.95)
  abline(v=c(-1,1), lty=2, col="grey50")
  abline(h=p05_fdr, lty=2, col="grey50")
  
  tops <- unique(c(head(tt[order(tt$adj.P.Val), "feature_id"], 10), qc_genes))
  tops <- tops[tops %in% tt$feature_id]
  sel  <- tt$feature_id %in% tops
  text(tt$logFC[sel], y[sel], labels = tt$feature_id[sel], pos=4, cex=0.8)
}

png(file.path(out_dir, "volcano.png"), width=1400, height=1100, res=150)
plot_with_margins(plot_volcano())
dev.off()

plot_ma <- function() {
  ttl <- wrap_title(sprintf("%s MA (Brueffer, n+=%d, n-=%d)", out_tag, n_pos, n_neg), 55)
  plot(fit$Amean, fit2$coefficients[,1],
       xlab="AveExpr", ylab="log2FC", main=ttl, pch=21, bg="white",
       col="black", cex.main = 0.95)
  abline(h=0, lty=2, col="grey50")
}

png(file.path(out_dir, "MA.png"), width=1400, height=1100, res=150)
plot_with_margins(plot_ma())
dev.off()

# ---------------- Concordance vs Wolf (FC2 vs FC2) ----------------
f_fc2 <- canon_fc2

if (file.exists(f_fc2)) {
  L <- fread(f_fc2, check.names = FALSE)
  stopifnot(all(c("feature_id", "logFC") %in% names(L)))
  L <- L[, .(feature_id, logFC)]
  
  wolf_dir <- file.path(dirname(script_dir), "wolf et al", "aim2_outputs", "GSE194040_ISPY2_mRNA_BPsubtype")
  wolf_opts <- c(
    file.path(wolf_dir, "FILT_HER2pos_vs_HER2neg_DE_FDR0.05_LFC1.tsv"),
    file.path(wolf_dir, "GSE194040_ISPY2_mRNA_BPsubtype_DE_filtered_FDRlt0.05_absLog2FCge1.tsv")
  )
  wolf_file <- wolf_opts[file.exists(wolf_opts)][1]
  
  if (length(wolf_file)) {
    W <- fread(wolf_file, check.names = FALSE)
    
    if ("gene_symbol" %in% names(W) && !"feature_id" %in% names(W)) W$feature_id <- W$gene_symbol
    if (!"logFC" %in% names(W)) {
      cand <- intersect(c("log2FC_HER2pos_vs_HER2neg", "log2FC", "LFC", "logFC"), names(W))
      if (length(cand)) setnames(W, cand[1], "logFC")
    }
    
    if (all(c("feature_id","logFC") %in% names(W))) {
      W_sub <- W[, .(feature_id, logFC_wolf = logFC)]
      M <- merge(L, W_sub, by = "feature_id")
      
      if (nrow(M)) {
        M[, same_direction := sign(logFC) == sign(logFC_wolf)]
        fc2_out  <- file.path(out_dir, "FILT_with_concordance_vs_Wolf_FC2.tsv")
        summ_out <- file.path(out_dir, "FILT_with_concordance_vs_Wolf_FC2_summary.tsv")
        fwrite(M, fc2_out, sep = "\t")
        
        summary_dt <- data.table(
          Brueffer_FC2_genes = nrow(L),
          Wolf_FC2_genes     = nrow(W_sub),
          Overlap_FC2_genes  = nrow(M),
          Same_direction     = sum(M$same_direction),
          Opposite_direction = nrow(M) - sum(M$same_direction),
          Pct_same_direction = 100 * sum(M$same_direction) / nrow(M)
        )
        fwrite(summary_dt, summ_out, sep = "\t")
        message("[Concordance] FC2 vs FC2 summary written: ", summ_out)
      } else {
        message("[Concordance] No overlap with Wolf.")
      }
    } else {
      message("[Concordance] Wolf file lacked feature_id/logFC after normalization.")
    }
  } else {
    message("[Concordance] Wolf FC2 file not found; skipping.")
  }
} else {
  message("[Concordance] Brueffer FC2 file not found; skipping.")
}

# ---------------- run log ----------------
erbb2_line <- { i <- which(tt$feature_id == "ERBB2")[1]
if (!is.na(i)) sprintf("ERBB2: log2FC=%.3f  FDR=%.3g", tt$logFC[i], tt$adj.P.Val[i]) else "ERBB2: NA"
}
grb7_line <- { i <- which(tt$feature_id == "GRB7")[1]
if (!is.na(i)) sprintf("GRB7:  log2FC=%.3f  FDR=%.3g", tt$logFC[i], tt$adj.P.Val[i]) else "GRB7: NA"
}

log_txt <- c(
  sprintf("RUN LOG – %s", out_tag),
  sprintf("Matrix: %s", basename(expr_file)),
  sprintf("Labels: %s", basename(lab_file)),
  sprintf("Samples kept: %d", length(grp)),
  sprintf("HER2 sizes: pos=%d neg=%d", n_pos, n_neg),
  sprintf("Wrote meta USED: %s", meta_used_tsv),
  sprintf("Contrast USED: %s", contrast_used_txt),
  sprintf("FDR<0.05: %d", nrow(tt_fdr)),
  sprintf("FDR<0.05 & |log2FC|>=1: %d", nrow(tabf)),
  "Outputs:",
  paste0("  - ", all_path),
  paste0("  - ", canon_fdr),
  paste0("  - ", canon_fc2),
  paste0("  - ", file.path(out_dir, "volcano.png")),
  paste0("  - ", file.path(out_dir, "MA.png")),
  paste0("  - ", file.path(out_dir, "CORR_HER2pos_vs_HER2neg_filtered.png")),
  paste0("  - ", file.path(out_dir, "FILT_with_concordance_vs_Wolf_FC2.tsv")),
  paste0("  - ", file.path(out_dir, "FILT_with_concordance_vs_Wolf_FC2_summary.tsv")),
  erbb2_line,
  grb7_line
)

writeLines(log_txt, file.path(out_dir, "RUN_LOG.txt"))
message("GSE81538 Brueffer RNA DE complete → ", out_dir)
# =======================================================================
