# her2_pipeline.R
# Core functions for running HER2+ vs HER2- differential expression (limma) across cohorts.
#
# This script:
#   - reads an expression matrix (features x samples)
#   - parses GEO series_matrix metadata to assign HER2 status
#   - maps feature IDs to gene symbols (biomaRt), keeping protein-coding + Swiss-Prot where possible
#   - runs limma and writes a standardized DE table + basic QC plots
#
# Note: biomaRt mapping requires internet access unless an ID map CSV is provided.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(dplyr)
  library(biomaRt)
})

# ----- small plotting helpers -----
wrap_title <- function(s, width = 60) paste(strwrap(s, width = width), collapse = "\n")
plot_with_margins <- function(expr) {
  op <- par(mar = c(5, 5, 5, 2))
  on.exit(par(op), add = TRUE)
  force(expr)
}

# ---------------- I/O + parsing helpers ----------------
read_expr_any <- function(path) {
  if (!file.exists(path)) stop("Missing expression file: ", path)
  dt <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  fid <- dt[[1]]
  fid[is.na(fid) | fid == ""] <- paste0("feat_", seq_len(sum(is.na(fid) | fid == "")))
  fid <- make.unique(as.character(fid))
  X <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(X) <- "double"
  rownames(X) <- fid
  X
}

collapse_reps <- function(E) {
  base <- function(x) {
    x <- gsub("[._-](rep|r)\\d+$", "", x, ignore.case = TRUE)
    x <- gsub("\\.\\d+$", "", x)
    x
  }
  bid <- base(colnames(E))
  idx <- split(seq_along(bid), bid)
  M <- vapply(
    idx,
    function(ix) {
      if (length(ix) == 1) E[, ix, drop = FALSE]
      else matrix(
        rowMeans(E[, ix, drop = FALSE], na.rm = TRUE),
        ncol = 1,
        dimnames = list(rownames(E), colnames(E)[ix][1])
      )
    },
    numeric(nrow(E))
  )
  out <- matrix(M, nrow = nrow(E))
  rownames(out) <- rownames(E)
  colnames(out) <- names(idx)
  out
}

# ---- series_matrix parsing ----
get_tablist <- function(series_path, token_regex) {
  hdr <- readLines(series_path, warn = FALSE)
  ln <- grep(token_regex, hdr, value = TRUE)[1]
  if (!length(ln) || is.na(ln)) return(NULL)
  parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  gsub('^"|"$', "", parts[-1])
}

get_all_chars <- function(series_path) {
  hdr <- readLines(series_path, warn = FALSE)
  grep("^!Sample_characteristics_ch1", hdr, value = TRUE)
}

label_from_raw <- function(x, pos = "Positive", neg = "Negative") {
  v <- tolower(trimws(x))
  out <- ifelse(grepl("amplified|amp\\b|\\b3\\+\\b|her2[ +]*\\+|positive", v), pos, NA_character_)
  out <- ifelse(grepl("non-?amplified|nonamp|\\b1\\+\\b|\\b0\\b|negative|her2[ -]*-", v), neg, out)
  out <- ifelse(grepl("equivocal|\\b2\\+\\b|not assessed|unknown|na", v), NA_character_, out)
  out
}

pick_best_her2_line <- function(series_path, pos = "Positive", neg = "Negative") {
  chr <- get_all_chars(series_path)
  her2_lines <- chr[grepl("her2|erbb2", tolower(chr))]
  if (!length(her2_lines)) stop("No HER2-related characteristics in: ", series_path)
  
  score_line <- function(line) {
    raw <- gsub('^"|"$', "", strsplit(line, "\t", fixed = TRUE)[[1]][-1])
    lab <- label_from_raw(raw, pos, neg)
    posn <- sum(lab == pos, na.rm = TRUE)
    negn <- sum(lab == neg, na.rm = TRUE)
    list(score = posn + negn + 2 * min(posn, negn), raw = raw, lab = lab)
  }
  S <- lapply(her2_lines, score_line)
  best_idx <- which.max(vapply(S, `[[`, numeric(1), "score"))
  best <- S[[best_idx]]
  list(
    raw = best$raw,
    lab = best$lab,
    token = sub("\t.*$", "", her2_lines[best_idx])
  )
}

get_gsms <- function(series_path) {
  hdr <- readLines(series_path, warn = FALSE)
  ln <- grep("^!Sample_geo_accession", hdr, value = TRUE)[1]
  if (!length(ln) || is.na(ln)) stop("No GSM line in: ", series_path)
  parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  gsub('^"|"$', "", parts[-1])
}

get_platforms <- function(series_path) {
  hdr <- readLines(series_path, warn = FALSE)
  ln <- grep("^!Sample_platform_id", hdr, value = TRUE)[1]
  if (!length(ln) || is.na(ln)) return(NULL)
  parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  gsub('^"|"$', "", parts[-1])
}

# ---------------- ID mapping (Swiss-Prot + protein-coding) ----------------
detect_id_type <- function(ids) {
  n_ens <- sum(grepl("^ENSG\\d+(\\.\\d+)?$", ids))
  n_sym <- sum(grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", ids))
  if (n_ens >= n_sym) "ensembl_gene_id" else "hgnc_symbol"
}

sanitize_ids <- function(ids) {
  ids <- as.character(ids)
  pipe_left <- sub("\\|.*$", "", ids)
  pipe_right <- sub("^.*\\|", "", ids)
  has_pipe <- grepl("\\|", ids)
  ids[has_pipe & grepl("^ENSG\\d+", pipe_right)] <- pipe_right[has_pipe & grepl("^ENSG\\d+", pipe_right)]
  ids[has_pipe & !grepl("^ENSG\\d+", pipe_right)] <- pipe_left[has_pipe & !grepl("^ENSG\\d+", pipe_right)]
  ids <- sub("\\.\\d+$", "", ids)
  trimws(ids)
}

build_id_map <- function(ids, cache_csv = NULL) {
  ids0 <- sanitize_ids(ids)
  id_type <- detect_id_type(ids0)
  mart <- suppressMessages(biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl"))
  attrs <- c(id_type, "hgnc_symbol", "uniprotswissprot", "gene_biotype")
  res <- biomaRt::getBM(attributes = attrs, filters = id_type, values = ids0, mart = mart)
  names(res) <- c("source_id", "gene_symbol", "uniprot_swissprot", "gene_biotype")
  res$protein_coding <- res$gene_biotype == "protein_coding"
  res <- dplyr::filter(res, !is.na(gene_symbol) & gene_symbol != "")
  if (!is.null(cache_csv)) utils::write.csv(res, cache_csv, row.names = FALSE)
  res
}

collapse_to_one_gene <- function(expr, map_df) {
  keep <- with(map_df, protein_coding & !is.na(uniprot_swissprot) & uniprot_swissprot != "")
  expr <- expr[keep, , drop = FALSE]
  map_df <- map_df[keep, , drop = FALSE]
  df <- data.frame(gene_symbol = map_df$gene_symbol, expr, check.names = FALSE)
  agg <- df |>
    dplyr::group_by(gene_symbol) |>
    dplyr::summarize(dplyr::across(where(is.numeric), \(x) median(x, na.rm = TRUE)), .groups = "drop")
  mat <- as.matrix(agg[, -1, drop = FALSE])
  rownames(mat) <- agg$gene_symbol
  list(expr = mat, map = map_df)
}

# ---------------- Main runner ----------------
her2_run <- function(dataset_dir,
                     expr_file,
                     series_files,
                     out_prefix,
                     collapse_replicates = FALSE,
                     pos_label = "Positive",
                     neg_label = "Negative",
                     her2_definition = "study-defined",
                     id_map_file = NULL) {
  
  message("Loading expression matrix")
  E <- read_expr_any(file.path(dataset_dir, expr_file))
  message("   first 5 feature IDs: ", paste(head(rownames(E), 5), collapse = ", "))
  
  # ---- detect whether log transform is needed ----
  transform_mode <- match.arg(getOption("her2.transform_mode", "auto"),
                              c("auto", "log2", "none"))
  vals <- as.numeric(E[is.finite(E)])
  p99 <- stats::quantile(vals, 0.99, na.rm = TRUE)
  mx <- max(vals, na.rm = TRUE)
  do_log <- switch(transform_mode,
                   "log2" = TRUE,
                   "none" = FALSE,
                   "auto" = { (p99 > 50) || (mx > 100) })
  if (do_log) {
    E <- log2(E + 1)
    message("   log2-transformed (mode=", transform_mode,
            "; p99=", round(p99, 2), ", max=", round(mx, 2), ")")
  } else {
    message("   assumed already on log scale (mode=", transform_mode,
            "; p99=", round(p99, 2), ", max=", round(mx, 2), ")")
  }
  
  if (collapse_replicates) {
    E <- collapse_reps(E)
    message("   collapsed technical replicates; samples=", ncol(E))
  }
  
  # ----- ID mapping & protein-coding Swiss-Prot filter -----
  message("ID mapping & protein-coding filter")
  row_ids <- rownames(E)
  if (!is.null(id_map_file) && file.exists(file.path(dataset_dir, id_map_file))) {
    map_df <- read.csv(file.path(dataset_dir, id_map_file), check.names = FALSE)
    stopifnot(all(c("source_id", "gene_symbol", "uniprot_swissprot", "protein_coding") %in% names(map_df)))
    map_df <- map_df[map_df$source_id %in% row_ids, , drop = FALSE]
  } else {
    cache_csv <- file.path(dataset_dir, paste0(out_prefix, "_id_map.csv"))
    map_df <- build_id_map(row_ids, cache_csv = cache_csv)
  }
  
  keep_rows <- intersect(row_ids, map_df$source_id)
  if (length(keep_rows) < 2000) warning("Only ", length(keep_rows), " rows mapped—check ID type.")
  E <- E[keep_rows, , drop = FALSE]
  
  # Deduplicate per source_id, preferring rows with Swiss-Prot where possible
  map_df <- map_df[map_df$source_id %in% rownames(E), , drop = FALSE]
  has_uni <- !is.na(map_df$uniprot_swissprot) & map_df$uniprot_swissprot != ""
  ord <- order(map_df$source_id, -as.integer(has_uni))
  map_df <- map_df[ord, , drop = FALSE]
  map_df <- map_df[!duplicated(map_df$source_id), , drop = FALSE]
  map_df <- map_df[match(rownames(E), map_df$source_id), , drop = FALSE]
  
  collapsed <- collapse_to_one_gene(E, map_df)
  E <- collapsed$expr
  map_df <- collapsed$map
  
  # ----- Parse series_matrix metadata and build a phenotype table keyed by GSM -----
  stopifnot(length(series_files) >= 1)
  
  pheno_list <- list()
  for (i in seq_along(series_files)) {
    spath <- file.path(dataset_dir, series_files[i])
    if (!file.exists(spath)) stop("Missing series matrix: ", spath)
    
    gsms <- get_gsms(spath)
    plats <- get_platforms(spath)
    
    # Make platform vector per sample
    if (is.null(plats)) {
      plat_vec <- rep(paste0("S", i), length(gsms))
    } else if (length(plats) == 1) {
      plat_vec <- rep(plats, length(gsms))
    } else if (length(plats) == length(gsms)) {
      plat_vec <- plats
    } else {
      major <- names(which.max(table(plats)))
      plat_vec <- rep(major, length(gsms))
    }
    
    titles <- get_tablist(spath, "^!Sample_title")
    if (is.null(titles)) titles <- rep(NA_character_, length(gsms))
    
    h <- pick_best_her2_line(spath, pos_label, neg_label)
    if (length(h$lab) != length(gsms)) stop("HER2 label length mismatch in: ", spath)
    
    pheno_list[[i]] <- data.frame(
      gsm = gsms,
      sample_title = titles,
      HER2_status = h$lab,
      HER2_raw = h$raw,
      platform = plat_vec,
      stringsAsFactors = FALSE
    )
    message("   ", basename(spath), "  HER2 field: ", h$token)
  }
  pheno <- do.call(rbind, pheno_list)
  
  # ----- output dirs -----
  prep_dir <- file.path(dataset_dir, paste0(out_prefix, "_prep"))
  out_dir <- file.path(dataset_dir, paste0("aim2_outputs/", out_prefix))
  dir.create(prep_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(pheno, file.path(prep_dir, "HER2_mapping_by_GSM.csv"), row.names = FALSE)
  
  # ----- Align expression columns to GEO metadata (GSM-first) -----
  colE <- trimws(colnames(E))
  pheno$gsm <- trimws(pheno$gsm)
  
  ov_gsm <- length(intersect(colE, pheno$gsm))
  message(sprintf("   overlap GSM: %d  (E=%d, pheno=%d)", ov_gsm, length(colE), nrow(pheno)))
  
  if (ov_gsm < 6) {
    # Fallback: align by title if expression columns are titles
    id_title <- intersect(colE, pheno$sample_title)
    if (length(id_title) >= 6) {
      message("   WARNING: expression columns match titles, not GSMs; aligning by title fallback")
      keep <- id_title
      E <- E[, keep, drop = FALSE]
      ann <- pheno[match(keep, pheno$sample_title), c("HER2_status", "platform"), drop = FALSE]
    } else {
      stop("Could not align expression to series metadata by GSM or title. Check expression column names.")
    }
  } else {
    keep <- intersect(colE, pheno$gsm)
    E <- E[, keep, drop = FALSE]
    ann <- pheno[match(keep, pheno$gsm), c("HER2_status", "platform"), drop = FALSE]
  }
  
  # Drop NA labels before limma
  na_n <- sum(is.na(ann$HER2_status))
  if (na_n > 0) message("   dropping ", na_n, " samples with NA/Equivocal HER2")
  keep2 <- !is.na(ann$HER2_status)
  E <- E[, keep2, drop = FALSE]
  ann <- ann[keep2, , drop = FALSE]
  
  # Group sizes + guardrails
  pos_n <- sum(ann$HER2_status == pos_label)
  neg_n <- sum(ann$HER2_status == neg_label)
  message(sprintf("   groups: HER2+ = %d, HER2- = %d", pos_n, neg_n))
  if (pos_n < 3 || neg_n < 3) stop("Groups too small after filtering (need ≥3 per arm).")
  
  # Save matched mapping
  write.csv(
    data.frame(sample = colnames(E),
               HER2_status = ann$HER2_status,
               platform = ann$platform),
    file.path(prep_dir, "HER2_mapping_matched.csv"),
    row.names = FALSE
  )
  
  # ----- QC: ERBB2/GRB7 boxplots -----
  qc_dir <- out_dir
  tgts <- c("ERBB2", "GRB7")
  present <- tgts[tgts %in% rownames(E)]
  n_pos_now <- sum(ann$HER2_status == pos_label, na.rm = TRUE)
  n_neg_now <- sum(ann$HER2_status == neg_label, na.rm = TRUE)
  
  if (length(present) > 0) {
    for (g in present) {
      fn_pdf <- file.path(qc_dir, sprintf("qc_%s_boxplot.pdf", g))
      fn_png <- file.path(qc_dir, sprintf("qc_%s_boxplot.png", g))
      ttl <- sprintf("%s %s | %s (n+=%d, n-=%d)", out_prefix, g, her2_definition, n_pos_now, n_neg_now)
      ttl <- wrap_title(ttl, 55)
      
      pdf(fn_pdf, width = 7.5, height = 6)
      plot_with_margins({
        boxplot(split(as.numeric(E[g, ]), ann$HER2_status),
                main = ttl, ylab = "Expression", xlab = "HER2 group",
                col = "grey90", border = "grey30", cex.main = 0.95)
      })
      dev.off()
      
      png(fn_png, width = 1400, height = 1100, res = 150)
      plot_with_margins({
        boxplot(split(as.numeric(E[g, ]), ann$HER2_status),
                main = ttl, ylab = "Expression", xlab = "HER2 group",
                col = "grey90", border = "grey30", cex.main = 0.95)
      })
      dev.off()
    }
  } else {
    warning("ERBB2/GRB7 not found after mapping—verify ID mapping and gene naming.")
  }
  
  # ----- limma -----
  zv <- apply(E, 1, sd, na.rm = TRUE) == 0
  if (any(zv)) E <- E[!zv, , drop = FALSE]
  
  grp <- factor(ifelse(ann$HER2_status == pos_label, "HER2pos", "HER2neg"),
                levels = c("HER2neg", "HER2pos"))
  platform <- factor(ann$platform)
  
  if (length(unique(platform)) > 1) {
    design <- model.matrix(~ 0 + grp + platform)
  } else {
    design <- model.matrix(~ 0 + grp)
  }
  colnames(design) <- make.names(colnames(design))
  colnames(design) <- sub("^grp", "", colnames(design))
  
  fit <- lmFit(E, design)
  ct <- makeContrasts(HER2pos_vs_HER2neg = HER2pos - HER2neg, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, ct), trend = TRUE, robust = TRUE)
  
  res <- topTable(fit2, coef = "HER2pos_vs_HER2neg", number = Inf, sort.by = "P")
  res$gene_symbol <- rownames(res)
  res$log2FC_HER2pos_vs_HER2neg <- res$logFC
  res$cohort_id <- out_prefix
  res$HER2_definition <- her2_definition
  res$n_pos <- sum(grp == "HER2pos")
  res$n_neg <- sum(grp == "HER2neg")
  
  res_out <- res[, c("gene_symbol", "log2FC_HER2pos_vs_HER2neg", "AveExpr", "t",
                     "P.Value", "adj.P.Val", "B", "cohort_id", "HER2_definition", "n_pos", "n_neg")]
  
  res_out_sig <- res_out[res_out$`adj.P.Val` < 0.05, , drop = FALSE]
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.table(res_out_sig,
              file.path(out_dir, paste0(out_prefix, "_DE_HER2pos_vs_HER2neg.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # ----- Volcano -----
  fnv_pdf <- file.path(out_dir, "volcano.pdf")
  fnv_png <- file.path(out_dir, "volcano.png")
  sig <- res$`adj.P.Val` < 0.05
  p05 <- -log10(0.05)
  
  plot_volcano <- function() {
    ttl <- paste0(out_prefix, " Volcano (", her2_definition, ", n+=",
                  sum(ann$HER2_status == pos_label), ", n-=",
                  sum(ann$HER2_status == neg_label), ")")
    ttl <- wrap_title(ttl, 55)
    
    with(res, plot(logFC, -log10(P.Value),
                   xlab = "log2FC (pos vs neg)", ylab = "-log10(P)",
                   main = ttl, pch = 21, bg = ifelse(sig, "black", "white"),
                   col = "black", cex.main = 0.95))
    abline(v = c(-1, 1), lty = 2, col = "grey50")
    abline(h = p05, lty = 2, col = "grey50")
    
    lab_genes <- unique(c(head(res[order(res$`adj.P.Val`), "gene_symbol"], 10), "ERBB2", "GRB7"))
    lab_genes <- lab_genes[lab_genes %in% res$gene_symbol]
    sel <- res$gene_symbol %in% lab_genes
    with(res[sel, ], text(logFC, -log10(P.Value), labels = gene_symbol, pos = 4, cex = 0.8))
  }
  pdf(fnv_pdf, width = 7.5, height = 6); plot_with_margins(plot_volcano()); dev.off()
  png(fnv_png, width = 1400, height = 1100, res = 150); plot_with_margins(plot_volcano()); dev.off()
  
  # ----- MA -----
  fnm_pdf <- file.path(out_dir, "MA.pdf")
  fnm_png <- file.path(out_dir, "MA.png")
  plot_ma <- function() {
    ttl <- paste0(out_prefix, " MA (", her2_definition, ", n+=",
                  sum(ann$HER2_status == pos_label), ", n-=",
                  sum(ann$HER2_status == neg_label), ")")
    ttl <- wrap_title(ttl, 55)
    plot(fit$Amean, fit2$coefficients[, "HER2pos_vs_HER2neg"],
         xlab = "AveExpr", ylab = "log2FC", main = ttl, pch = 21,
         bg = "white", col = "black", cex.main = 0.95)
    abline(h = 0, lty = 2, col = "grey50")
  }
  pdf(fnm_pdf, width = 7.5, height = 6); plot_with_margins(plot_ma()); dev.off()
  png(fnm_png, width = 1400, height = 1100, res = 150); plot_with_margins(plot_ma()); dev.off()
  
  invisible(list(E = E, ann = ann, design = design, meta = pheno))
}

write_run_log <- function(base, out_prefix) {
  out_dir <- file.path(base, "aim2_outputs", out_prefix)
  out_tsv <- file.path(out_dir, paste0(out_prefix, "_DE_HER2pos_vs_HER2neg.tsv"))
  if (!file.exists(out_tsv)) stop("DE table not found: ", out_tsv)
  
  x <- read.delim(out_tsv, check.names = FALSE)
  
  pick_line <- function(df, g) {
    row <- df[df$gene_symbol == g, , drop = FALSE]
    if (nrow(row)) {
      sprintf("%s: log2FC=%.3f, FDR=%.3g", g,
              row$log2FC_HER2pos_vs_HER2neg[1], row$`adj.P.Val`[1])
    } else sprintf("%s: NA", g)
  }
  
  erbb2_line <- pick_line(x, "ERBB2")
  grb7_line <- pick_line(x, "GRB7")
  
  n_sig <- nrow(x)
  n_pos <- unique(x$n_pos)
  n_neg <- unique(x$n_neg)
  her2_def <- unique(x$HER2_definition)
  
  transform_mode <- getOption("her2.transform_mode", default = "auto")
  
  log_txt <- c(
    sprintf("RUN LOG — %s", out_prefix),
    sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    sprintf("Cohort ID: %s", out_prefix),
    sprintf("HER2 definition: %s", her2_def),
    sprintf("Group sizes: n+=%s, n-=%s", paste(n_pos, collapse = ","), paste(n_neg, collapse = ",")),
    "",
    sprintf("Expression transform mode: %s", transform_mode),
    "Mapping: Ensembl/symbol → gene_symbol via biomaRt; filtered to protein-coding + Swiss-Prot.",
    "Deduplication: prefer rows with Swiss-Prot per source_id.",
    "",
    sprintf("FDR < 0.05 hits: %d", n_sig),
    sprintf("QC targets: %s; %s", erbb2_line, grb7_line),
    "",
    "Key outputs:",
    sprintf("- DE table: %s", out_tsv),
    sprintf("- Volcano: %s", file.path(out_dir, "volcano.pdf")),
    sprintf("- MA: %s", file.path(out_dir, "MA.pdf")),
    sprintf("- QC ERBB2 boxplot: %s",
            if (file.exists(file.path(out_dir, "qc_ERBB2_boxplot.pdf")))
              file.path(out_dir, "qc_ERBB2_boxplot.pdf") else "missing"),
    sprintf("- QC GRB7 boxplot: %s",
            if (file.exists(file.path(out_dir, "qc_GRB7_boxplot.pdf")))
              file.path(out_dir, "qc_GRB7_boxplot.pdf") else "missing"),
    "",
    "Session:",
    capture.output(sessionInfo())
  )
  
  fn <- file.path(out_dir, "RUN_LOG.txt")
  writeLines(log_txt, fn)
  message("Wrote: ", fn)
  invisible(fn)
}

# ======================== end her2_pipeline.R ========================
