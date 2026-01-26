# ============================ wolf_harps_pipeline.R ============================
# Cohort: HArPS labeling (RPS-5), RNA DE (HArPS+ vs HArPS-), optional concordance
# against HER2-status DE, and optional RPPA/proteome DE blocks.
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

# Core inputs
SUPP3_XLSX <- file.path(script_dir, "NIHMS1829047-supplement-3.xlsx")
RNA_FILE   <- file.path(script_dir, "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt")

# Optional molecular layers (only used if toggled on and files exist)
RPPA_FILE  <- file.path(script_dir, "GSE196093_RPPA_matrix.tsv")
PROT_FILE  <- file.path(script_dir, "wolf_proteome_expr.tsv")

# Output base
OUT_ROOT <- file.path(script_dir, "aim2_outputs")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# Force all RNA outputs here
RNA_OUT_DIR <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/Wolf_HARP_RNA"
dir.create(RNA_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Thresholds (used for filtered significant sets)
FDR_MAX <- 0.05
LFC_MIN <- 1.0

# Toggles
RUN_RNA        <- TRUE
RUN_CONCORD    <- TRUE   # requires HER2_DE_FILE
MAKE_PLOTS     <- TRUE
RUN_RPPA       <- FALSE
RUN_PROTEOME   <- FALSE

# HER2-status DE input (produced elsewhere)
HER2_DE_FILE <- file.path(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al",
  "aim2_outputs",
  "GSE194040_ISPY2_mRNA_BPsubtype",
  "ALL_HER2pos_vs_HER2neg_DE.tsv"
)

# Concordance output folder
CONCORD_OUT_DIR <- RNA_OUT_DIR

# Optional prioritized list file (used for overlap subset plots/tables)
P76_FILE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv"

# ==============================================================================
# HELPERS: labeling + IO
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

# ---------- wide reader (repairs header drift) ----------
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
  
  # coerce numeric where possible
  for (nm in samp_cols) {
    if (!is.numeric(dt[[nm]])) {
      v <- as.character(dt[[nm]]); v[is.na(v)] <- ""
      suppressWarnings(num <- as.numeric(v))
      if (sum(nchar(v) > 0 & is.na(num)) == 0) dt[[nm]] <- num
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
  E
}

# ---------- LIMMA ----------
run_limma <- function(expr, meta_dt, group_col, pos_label = "HARP+", neg_label = "HARP-") {
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
    feature    = rownames(tt),
    feature_id = rownames(tt),
    logFC      = tt$logFC,
    AveExpr    = tt$AveExpr,
    t          = tt$t,
    P.Value    = tt$P.Value,
    adj.P.Val  = tt$adj.P.Val
  )
}

filt_hits <- function(dt, fdr_max, lfc_min) {
  dt[!is.na(adj.P.Val) & adj.P.Val < fdr_max & abs(logFC) >= lfc_min]
}

# ==============================================================================
# HELPERS: correlation + scatter styling (open circles + lm line + annotation)
# ==============================================================================

spearman_stats <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value)
}

fmt_p <- function(p) {
  if (is.na(p)) return("=NA")
  if (p < 1e-4) return("<1E-4")
  paste0("=", sprintf("%.1E", p))
}

make_scatter_bateman_style <- function(dt, xcol, ycol, title, out_png,
                                       xlab = "her2_logFC", ylab = "harp_logFC",
                                       w = 7.2, h = 5.0, dpi = 300,
                                       title_size = 18) {
  stopifnot(nrow(dt) > 2)
  
  st  <- spearman_stats(dt[[xcol]], dt[[ycol]])
  ann <- sprintf("Rho=%.3f, P%s", st$rho, fmt_p(st$p))
  
  p <- ggplot(dt, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(shape = 1, size = 2.2, alpha = 0.9) +      # open circles
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    theme_classic(base_size = 16) +
    labs(title = title, x = xlab, y = ylab) +
    theme(
      plot.title = element_text(size = title_size, hjust = 0.5),
      plot.margin = margin(12, 18, 12, 12)
    ) +
    annotate("text",
             x = Inf, y = -Inf,
             label = ann,
             hjust = 1.05, vjust = -0.6,
             size = 6)
  
  ggsave(out_png, p, width = w, height = h, dpi = dpi)
  list(plot = p, rho = st$rho, p = st$p)
}

# ==============================================================================
# STEP 1: metadata (RPS-5 → HER2 + HArPS)
# ==============================================================================

stopifnot(file.exists(SUPP3_XLSX))
supp3 <- as.data.table(read_xlsx(SUPP3_XLSX, sheet = 1))
meta  <- supp3[, .(`Patient Identifier`, `RPS-5`)]
setnames(meta, c("Patient Identifier","RPS-5"), c("sample_id","RPS5"))

meta[, HER2_status := her2_status_from_rps5(RPS5)]
meta[, HARP        := map_harp_from_rps5(RPS5)]

meta_labels <- meta[, .(sample_id, HER2_status, HARP, RPS5)]
fwrite(meta_labels, file.path(OUT_ROOT, "labels_from_RPS5.csv"))

meta_harp <- meta[!is.na(HARP), .(sample_id, HARP, RPS5)]
fwrite(meta_harp, file.path(OUT_ROOT, "metadata_HARP.csv"))

cat("HArPS group counts:\n")
print(table(meta_harp$HARP, useNA = "ifany"))

# ==============================================================================
# STEP 2: composition plot (optional)
# ==============================================================================

if (MAKE_PLOTS) {
  plot_dt <- copy(meta_labels)[!is.na(HER2_status)]
  plot_dt[, HARP_plot := fifelse(is.na(HARP), "Not assigned (HER2+)", as.character(HARP))]
  plot_dt[, HARP_plot := factor(HARP_plot, levels = c("HARP-","HARP+","Not assigned (HER2+)"))]
  
  p_comp <- ggplot(plot_dt, aes(x = HER2_status, fill = HARP_plot)) +
    geom_bar(position = "fill", width = 0.85) +
    theme_classic(base_size = 12) +
    labs(
      title = "Cohort composition by HER2 status and HArPS label",
      x = "HER2 status",
      y = "Proportion of samples",
      fill = "HArPS label"
    )
  
  ggsave(
    filename = file.path(OUT_ROOT, "composition_by_HER2_status.png"),
    plot     = p_comp,
    width    = 7, height = 5, dpi = 300
  )
}

# ==============================================================================
# STEP 3: RNA DE (HArPS+ vs HArPS-)
# ==============================================================================

E_rna <- NULL
rna_tt <- NULL

if (RUN_RNA) {
  stopifnot(file.exists(RNA_FILE))
  cat("\n--- Running LIMMA on RNA (HArPS+ vs HArPS-) ---\n")
  
  E_rna  <- read_expr_wide(RNA_FILE, feature_name = "gene_symbol")
  rna_tt <- run_limma(E_rna, meta_harp, "HARP")
  
  fwrite(rna_tt, file.path(RNA_OUT_DIR, "LIMMA_all.csv"))
  fwrite(rna_tt[adj.P.Val < FDR_MAX], file.path(RNA_OUT_DIR, "LIMMA_FDR0.05.csv"))
  fwrite(filt_hits(rna_tt, FDR_MAX, LFC_MIN), file.path(RNA_OUT_DIR, "LIMMA_FDR0.05_LFC1.csv"))
  
  cat("RNA DE complete:\n")
  cat("  all =", nrow(rna_tt), "\n")
  cat("  FDR<0.05 =", nrow(rna_tt[adj.P.Val < FDR_MAX]), "\n")
  cat("  FDR<0.05 & |logFC|>=1 =", nrow(filt_hits(rna_tt, FDR_MAX, LFC_MIN)), "\n")
}

# ==============================================================================
# STEP 4: Concordance vs HER2-status DE (relationship plot)
# ==============================================================================

if (RUN_CONCORD) {
  
  if (!file.exists(HER2_DE_FILE)) {
    warning("HER2-status DE file not found; skipping:\n  ", HER2_DE_FILE)
  } else if (is.null(rna_tt) || !nrow(rna_tt)) {
    warning("HArPS RNA DE results not available; skipping concordance.")
  } else {
    
    her2 <- fread(HER2_DE_FILE)
    harp <- copy(rna_tt)
    
    normalize_feature_col <- function(dt, candidates = c("feature_id","feature","gene","gene_symbol","GeneSymbol","symbol")) {
      hit <- candidates[candidates %in% names(dt)][1]
      if (is.na(hit)) hit <- names(dt)[1]
      setnames(dt, hit, "feature_id")
      dt
    }
    
    her2 <- normalize_feature_col(her2)
    harp <- normalize_feature_col(harp, candidates = c("feature_id","feature","gene","gene_symbol","symbol"))
    
    req_cols <- c("feature_id","logFC","P.Value","adj.P.Val")
    if (!all(req_cols %in% names(her2))) stop("HER2 DE file missing: ", paste(setdiff(req_cols, names(her2)), collapse = ", "))
    if (!all(req_cols %in% names(harp))) stop("HArPS DE table missing: ", paste(setdiff(req_cols, names(harp)), collapse = ", "))
    
    her2_sub <- her2[, .(feature_id, her2_logFC = logFC, her2_P = P.Value, her2_FDR = adj.P.Val)]
    harp_sub <- harp[, .(feature_id, harp_logFC = logFC, harp_P = P.Value, harp_FDR = adj.P.Val)]
    
    concord <- merge(her2_sub, harp_sub, by = "feature_id", all = FALSE)
    concord_dt <- as.data.table(concord)
    
    # Flags
    concord_dt[, both_FDR := (her2_FDR < FDR_MAX & harp_FDR < FDR_MAX)]
    concord_dt[, same_direction := (sign(her2_logFC) == sign(harp_logFC) & her2_logFC != 0 & harp_logFC != 0)]
    
    # Save full merged table
    fwrite(concord_dt, file.path(CONCORD_OUT_DIR, "HER2_vs_HARP_concordance_ALL.csv"))
    
    # Sanity checks to catch sign/contrast flips quickly
    cat("\nConcordance checks:\n")
    cat("  merged genes =", nrow(concord_dt), "\n")
    cat("  FDR<0.05 in both =", sum(concord_dt$both_FDR), "\n")
    cat("  same-direction among FDR<0.05 both =",
        sum(concord_dt$both_FDR & concord_dt$same_direction), "\n")
    cat("  opposite-direction among FDR<0.05 both =",
        sum(concord_dt$both_FDR & !concord_dt$same_direction), "\n")
    
    # --- Plot A: relationship using genes significant in BOTH (includes opposite-direction) ---
    plot_bothFDR <- concord_dt[both_FDR == TRUE]
    if (MAKE_PLOTS && nrow(plot_bothFDR) > 2) {
      make_scatter_bateman_style(
        dt = plot_bothFDR,
        xcol = "her2_logFC",
        ycol = "harp_logFC",
        title = "RNA: relationship (FDR<0.05 in both)",
        out_png = file.path(CONCORD_OUT_DIR, "HER2_vs_HARP_relationship_FDR0.05_in_both.png"),
        xlab = "her2_logFC",
        ylab = "harp_logFC",
        w = 7.6, h = 5.2, dpi = 300,
        title_size = 18
      )
    }
    
    # --- Plot B: relationship using genes significant in BOTH + same direction (often matches higher rho) ---
    plot_bothFDR_same <- concord_dt[both_FDR == TRUE & same_direction == TRUE]
    fwrite(plot_bothFDR_same, file.path(CONCORD_OUT_DIR, "HER2_vs_HARP_FDR0.05_in_both_same_direction.csv"))
    
    if (MAKE_PLOTS && nrow(plot_bothFDR_same) > 2) {
      make_scatter_bateman_style(
        dt = plot_bothFDR_same,
        xcol = "her2_logFC",
        ycol = "harp_logFC",
        title = "RNA: relationship (FDR<0.05 in both; same direction)",
        out_png = file.path(CONCORD_OUT_DIR, "HER2_vs_HARP_relationship_FDR0.05_in_both_same_direction.png"),
        xlab = "her2_logFC",
        ylab = "harp_logFC",
        w = 7.6, h = 5.2, dpi = 300,
        title_size = 16
      )
    }
    
    # --- 3-way overlap with prioritized list (if available) ---
    if (file.exists(P76_FILE)) {
      p76 <- fread(P76_FILE)
      possible_gene_cols <- c("gene_symbol","GeneSymbol","symbol","gene","feature_id","Feature","Gene")
      gene_col <- possible_gene_cols[possible_gene_cols %in% names(p76)][1]
      
      if (is.na(gene_col)) {
        warning("No gene column found in prioritized list. Columns: ", paste(names(p76), collapse = ", "))
      } else {
        p76_genes <- unique(as.character(p76[[gene_col]]))
        
        overlap_3way <- plot_bothFDR_same[feature_id %in% p76_genes]
        fwrite(overlap_3way, file.path(CONCORD_OUT_DIR, "Overlap_between_Wolf_HARP_HER2_and_P76.csv"))
        
        overlap_names <- overlap_3way[, .(gene = feature_id)]
        fwrite(overlap_names, file.path(CONCORD_OUT_DIR, "Overlap_between_Wolf_HARP_HER2_and_P76_geneList.csv"))
        
        # Optional: quick printed list for copy/paste into a slide
        cat("\nOverlap between Wolf HARP, HER2 & P76 (N=", nrow(overlap_names), "):\n", sep = "")
        if (nrow(overlap_names) > 0) print(overlap_names)
      }
    } else {
      message("Prioritized list file not found; skipping 3-way overlap:\n  ", P76_FILE)
    }
  }
}

# ==============================================================================
# STEP 5: Optional RPPA / Proteome DE (same HArPS labels)
# ==============================================================================

if (RUN_RPPA) {
  cat("\n--- Attempting LIMMA on RPPA (optional) ---\n")
  if (!file.exists(RPPA_FILE)) {
    message("  RPPA file not found: ", RPPA_FILE)
  } else {
    E_rppa <- tryCatch(read_expr_wide(RPPA_FILE, feature_name = "feature"),
                       error = function(e) { message("  Skipping RPPA: ", e$message); NULL })
    if (!is.null(E_rppa)) {
      rppa_tt <- run_limma(E_rppa, meta_harp, "HARP")
      rppa_dir <- file.path(OUT_ROOT, "HARP_RPPA")
      dir.create(rppa_dir, showWarnings = FALSE, recursive = TRUE)
      
      fwrite(rppa_tt, file.path(rppa_dir, "LIMMA_all.csv"))
      fwrite(rppa_tt[adj.P.Val < FDR_MAX], file.path(rppa_dir, "LIMMA_FDR0.05.csv"))
      fwrite(filt_hits(rppa_tt, FDR_MAX, LFC_MIN), file.path(rppa_dir, "LIMMA_FDR0.05_LFC1.csv"))
      cat("RPPA DE complete -> ", rppa_dir, "\n")
    }
  }
}

if (RUN_PROTEOME) {
  cat("\n--- Attempting LIMMA on proteome (optional) ---\n")
  if (!file.exists(PROT_FILE)) {
    message("  Proteome file not found: ", PROT_FILE)
  } else {
    E_prot <- tryCatch(read_expr_wide(PROT_FILE, feature_name = "feature"),
                       error = function(e) { message("  Skipping proteome: ", e$message); NULL })
    if (!is.null(E_prot)) {
      prot_tt <- run_limma(E_prot, meta_harp, "HARP")
      prot_dir <- file.path(OUT_ROOT, "HARP_Proteome")
      dir.create(prot_dir, showWarnings = FALSE, recursive = TRUE)
      
      fwrite(prot_tt, file.path(prot_dir, "LIMMA_all.csv"))
      fwrite(prot_tt[adj.P.Val < FDR_MAX], file.path(prot_dir, "LIMMA_FDR0.05.csv"))
      fwrite(filt_hits(prot_tt, FDR_MAX, LFC_MIN), file.path(prot_dir, "LIMMA_FDR0.05_LFC1.csv"))
      cat("Proteome DE complete -> ", prot_dir, "\n")
    }
  }
}

cat("\n✅ Finished. Outputs -> ", RNA_OUT_DIR, "\n")
# ==============================================================================
