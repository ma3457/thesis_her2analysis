# ========== run_RK_PROT.R — GLOBAL proteome DEP (HER2+/−) — FINAL (PNG-ONLY + -log10(FDR) + META USED/CONTRAST + CONCORDANCE FIXED) ==========
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(limma)
  library(stringr)
})

# -------- CONFIG --------
base_dir  <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/raj kumar et al"
prot_file <- file.path(base_dir, "proteomics_data", "20220330_APOLLO4C_GlobalProteomics.txt")
clin_xlsx <- file.path(base_dir, "APOLLO4C_ClinicalData_Freeze_V3_20220316.xlsx")
cohort_id <- "RajKumar_PROT_GLOBAL"

map_csv   <- file.path(base_dir, "Processed_filtered_matrix", "RK_map_PROT_to_APOLLO.csv")
out_dir   <- file.path(base_dir, "aim2_outputs", cohort_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_per_group <- 3

# Meta outputs (USED + contrast)  [template optional but kept; won't break if present]
meta_template_csv <- file.path(out_dir, paste0(cohort_id, "_meta_template.csv"))
meta_used_tsv     <- file.path(out_dir, paste0(cohort_id, "_meta_USED.tsv"))
contrast_used_txt <- file.path(out_dir, paste0(cohort_id, "_contrast_USED.txt"))

# Where to search for transcript DE outputs to mark concordance
search_roots <- c(
  base_dir,
  dirname(base_dir),
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA"
)

`%||%` <- function(a,b) if (!is.null(a)) a else b

# -------- HELPERS --------
normalize_id <- function(v){
  v <- trimws(as.character(v))
  v <- gsub("\\s+", "", v)
  v <- gsub("\\.", "-", v)   # AP.B3K8.ML -> AP-B3K8-ML
  toupper(v)
}

read_matrix_any <- function(path){
  x <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  ann_cols  <- intersect(names(x), c("Protein","Gene","Accession","Description"))
  samp_cols <- setdiff(names(x), ann_cols)
  
  for (nm in samp_cols) if (!is.numeric(x[[nm]])) x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
  
  feat <- if ("Gene" %in% names(x)) x$Gene else if ("Protein" %in% names(x)) x$Protein else x[[1]]
  feat <- trimws(as.character(feat))
  
  if (!("Gene" %in% names(x)) && "Protein" %in% names(x)) {
    guess <- sub("^.*\\|", "", x$Protein)
    guess <- sub("_HUMAN.*$", "", guess, ignore.case = TRUE)
    idx <- which(is.na(feat) | feat == "")
    feat[idx] <- guess[idx]
  }
  
  feat <- toupper(gsub("[^A-Za-z0-9_.-]", "", feat))
  M <- as.matrix(x[, samp_cols, drop = FALSE]); storage.mode(M) <- "numeric"
  rownames(M) <- make.unique(feat)
  
  keep <- !is.na(rownames(M)) & rownames(M) != ""
  M[keep, , drop = FALSE]
}

log2_auto <- function(M){
  vals <- as.numeric(M[is.finite(M)])
  if (!length(vals)) return(M)
  p99  <- stats::quantile(vals, 0.99, na.rm = TRUE)
  mx   <- suppressWarnings(max(vals, na.rm = TRUE))
  if (is.finite(p99) && (p99 > 50 || mx > 100)) { message(">>> log2-transforming proteome"); log2(M + 1) } else M
}

label_HER2 <- function(ihc, ratio){
  ihc_l <- trimws(toupper(gsub("\\s+","",as.character(ihc))))
  ratio <- suppressWarnings(as.numeric(ratio))
  ifelse(ihc_l=="3+","HER2pos",
         ifelse(ihc_l %in% c("0","1+"),"HER2neg",
                ifelse(ihc_l=="2+", ifelse(is.finite(ratio)&ratio>=2,"HER2pos",
                                           ifelse(is.finite(ratio)&ratio<2,"HER2neg",NA)),
                       ifelse(is.finite(ratio)&ratio>=2,"HER2pos",
                              ifelse(is.finite(ratio)&ratio<2,"HER2neg",NA)))))
}

robust_impute_and_filter <- function(M){
  all_bad <- rowSums(!is.finite(M)) == ncol(M)
  if (any(all_bad)) M <- M[!all_bad, , drop = FALSE]
  M[!is.finite(M)] <- NA
  
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

# -------- LOAD DATA --------
M <- read_matrix_any(prot_file)
M <- log2_auto(M)

clin_raw <- as.data.frame(read_xlsx(clin_xlsx))

apollo_col <- names(clin_raw)[grepl("APOLLO", names(clin_raw), ignore.case = TRUE) &
                                grepl("Tumor",  names(clin_raw), ignore.case = TRUE)]
if (!length(apollo_col)) apollo_col <- names(clin_raw)[grepl("APOLLO", names(clin_raw), ignore.case = TRUE)]
stopifnot(length(apollo_col) >= 1)

clin <- data.frame(
  APOLLO = normalize_id(clin_raw[[apollo_col[1]]]),
  HER2   = label_HER2(clin_raw$`HER2 IHC RESULT` %||% clin_raw$HER2.IHC.RESULT,
                      clin_raw$`HER2 RATIO`      %||% clin_raw$HER2.RATIO),
  stringsAsFactors = FALSE
)

# Template-style meta table derived from clinical sheet (safe; optional)
clin_meta <- data.table(
  SampleID_raw = normalize_id(clin_raw[[apollo_col[1]]]),
  SampleID     = normalize_id(clin_raw[[apollo_col[1]]]),
  HER2_IHC     = as.character(clin_raw$`HER2 IHC RESULT` %||% clin_raw$HER2.IHC.RESULT),
  HER2_RATIO   = as.character(clin_raw$`HER2 RATIO`      %||% clin_raw$HER2.RATIO),
  HER2_group   = label_HER2(clin_raw$`HER2 IHC RESULT` %||% clin_raw$HER2.IHC.RESULT,
                            clin_raw$`HER2 RATIO`      %||% clin_raw$HER2.RATIO),
  stringsAsFactors = FALSE
)

maybe_add <- function(colname, newname = colname) {
  if (colname %in% names(clin_raw)) clin_meta[[newname]] <<- as.character(clin_raw[[colname]])
}
maybe_add("PATIENT_ID", "PatientID")
maybe_add("SAMPLE_ID",  "SampleID_alt")
maybe_add("TUMOR_TYPE", "TumorType")

if (!file.exists(meta_template_csv)) {
  tmp <- copy(clin_meta)
  tmp[, Notes := ""]
  fwrite(tmp, meta_template_csv)
  message("Wrote meta template: ", meta_template_csv)
} else {
  message("Meta template exists: ", meta_template_csv)
}

# -------- ALIGNMENT --------
colnames(M) <- normalize_id(colnames(M))
keep <- intersect(colnames(M), clin$APOLLO)

if (length(keep) == 0) {
  if (!file.exists(map_csv)) {
    tmpl <- data.frame(prot_colname = colnames(M), APOLLO_ID = "", stringsAsFactors = FALSE)
    data.table::fwrite(tmpl, map_csv)
    stop(paste0("No overlap between proteome columns and clinical IDs. Wrote mapping template:\n", map_csv,
                "\nFill APOLLO_ID for each proteome column, then re-run."))
  }
  mp <- data.table::fread(map_csv, data.table = FALSE)
  stopifnot(all(c("prot_colname","APOLLO_ID") %in% names(mp)))
  mp$prot_colname <- normalize_id(mp$prot_colname)
  mp$APOLLO_ID    <- normalize_id(mp$APOLLO_ID)
  ren <- setNames(mp$APOLLO_ID, mp$prot_colname)
  colnames(M) <- ifelse(colnames(M) %in% names(ren) & ren[colnames(M)] != "", ren[colnames(M)], colnames(M))
  keep <- intersect(colnames(M), clin$APOLLO)
  if (length(keep) == 0) stop("Still 0 overlap after applying RK_map_PROT_to_APOLLO.csv")
}

M   <- M[, keep, drop = FALSE]
grp <- factor(clin$HER2[match(colnames(M), clin$APOLLO)], levels = c("HER2neg","HER2pos"))
ok  <- !is.na(grp)
M   <- M[, ok, drop = FALSE]
grp <- droplevels(grp[ok])

if (length(levels(grp)) < 2) stop("Grouping has <2 levels after filtering.")
if (min(table(grp)) < min_per_group) stop("Groups too small: ", paste(table(grp), collapse=" / "))

# ---- robust NA/variance handling
M <- robust_impute_and_filter(M)
stopifnot(nrow(M) >= 2)

# -------- EXPORT META USED + CONTRAST USED --------
meta_used <- data.table(
  SampleID = colnames(M),
  Group    = as.character(grp)
)

tmp_clin <- copy(clin_meta)
tmp_clin[, SampleID := normalize_id(SampleID)]
meta_used <- merge(meta_used, tmp_clin, by = "SampleID", all.x = TRUE)

if (file.exists(map_csv)) {
  mp <- fread(map_csv)
  if (all(c("prot_colname","APOLLO_ID") %in% names(mp))) {
    mp[, prot_colname := normalize_id(prot_colname)]
    mp[, APOLLO_ID    := normalize_id(APOLLO_ID)]
    setnames(mp, c("prot_colname","APOLLO_ID"), c("ProteomeColumn","SampleID"))
    meta_used <- merge(meta_used, mp, by = "SampleID", all.x = TRUE)
  }
}

fwrite(meta_used, meta_used_tsv, sep = "\t", quote = FALSE, na = "NA")

writeLines(c(
  "Mode: HER2",
  "Contrast: HER2pos - HER2neg",
  paste0("Cohort: ", cohort_id),
  paste0("n_pos: ", sum(grp == "HER2pos")),
  paste0("n_neg: ", sum(grp == "HER2neg"))
), contrast_used_txt)

message("Exported: ", meta_used_tsv)
message("Exported: ", contrast_used_txt)

# -------- LIMMA DEP --------
design <- model.matrix(~ 0 + grp)
colnames(design) <- c("HER2neg","HER2pos")

fit  <- lmFit(M, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(HER2pos - HER2neg, levels = design)),
               trend = TRUE, robust = TRUE)

tt <- topTable(fit2, number = Inf, sort.by = "P")
tt$feature_id <- rownames(tt)

# -------- OUTPUTS (full + canonical filtered: FDR, FDR+FC2) --------
full_cols <- c("feature_id","logFC","AveExpr","t","P.Value","adj.P.Val","B")
tt_full <- as.data.frame(tt)[, full_cols, drop = FALSE]
fwrite(tt_full, file.path(out_dir, "DE_full.tsv"), sep = "\t")

out_base <- paste0(cohort_id, "_DE_HER2pos_vs_HER2neg")
f_fdr <- file.path(out_dir, paste0(out_base, "_DE_filtered_FDRlt0.05.tsv"))
f_fc2 <- file.path(out_dir, paste0(out_base, "_DE_filtered_FDRlt0.05_absLog2FCge1.tsv"))

sig_FDR  <- tt[tt$adj.P.Val < 0.05, ]
sig_FDR2 <- sig_FDR[abs(sig_FDR$logFC) >= 1, ]

mk_out <- function(d){
  if (!nrow(d)) return(d)
  data.frame(
    feature_id = d$feature_id,
    log2FC_HER2pos_vs_HER2neg = d$logFC,
    AveExpr = d$AveExpr, t = d$t,
    P.Value = d$P.Value, adj.P.Val = d$adj.P.Val, B = d$B,
    cohort_id = cohort_id,
    HER2_definition = "IHC primary; ratio fallback (>=2.0)",
    n_pos = sum(grp == "HER2pos"),
    n_neg = sum(grp == "HER2neg"),
    check.names = FALSE
  )
}

fwrite(mk_out(sig_FDR),  f_fdr, sep = "\t")
fwrite(mk_out(sig_FDR2), f_fc2, sep = "\t")

# -------- ERBB2 / GRB7 QC (PNG-only) --------
qc_genes <- c("ERBB2", "GRB7")
qc_found <- intersect(qc_genes, rownames(M))

qc_rows <- list()
for (g in qc_genes) {
  if (!(g %in% qc_found)) { message("QC note: ", g, " not in rownames."); next }
  
  vals <- as.numeric(M[g, ])
  v_pos <- vals[grp == "HER2pos"]; v_neg <- vals[grp == "HER2neg"]
  m_pos <- mean(v_pos, na.rm = TRUE); s_pos <- sd(v_pos, na.rm = TRUE)
  m_neg <- mean(v_neg, na.rm = TRUE); s_neg <- sd(v_neg, na.rm = TRUE)
  log2fc <- m_pos - m_neg
  ttst <- try(stats::t.test(v_pos, v_neg, var.equal = FALSE), silent = TRUE)
  pval <- if (inherits(ttst, "try-error")) NA_real_ else ttst$p.value
  
  png(file.path(out_dir, paste0("qc_", g, "_boxplot.png")), width = 1200, height = 1000, res = 150)
  par(mar = c(5,5,4.5,2))
  boxplot(split(vals, grp),
          main = sprintf("%s (global proteome, log-scale)", g),
          ylab = "Expression", xlab = "HER2 group",
          col = c("grey85","grey85"), border = "grey25")
  stripchart(split(vals, grp), vertical = TRUE, method = "jitter",
             pch = 21, bg = "white", col = "black", add = TRUE)
  mtext(sprintf("mean±SD (pos)=%.2f±%.2f | (neg)=%.2f±%.2f | log2FC=%.2f | P=%.3g",
                m_pos, s_pos, m_neg, s_neg, log2fc, pval),
        side = 3, line = 0.5, cex = 0.85)
  dev.off()
  
  qc_rows[[g]] <- data.frame(
    gene = g,
    n_pos = length(v_pos), mean_pos = m_pos, sd_pos = s_pos,
    n_neg = length(v_neg), mean_neg = m_neg, sd_neg = s_neg,
    log2FC_pos_vs_neg = log2fc,
    ttest_pvalue = pval,
    stringsAsFactors = FALSE
  )
}

if (length(qc_rows)) {
  qc_tab <- do.call(rbind, qc_rows)
  fwrite(qc_tab, file = file.path(out_dir, "QC_ERBB2_GRB7_summary.tsv"), sep = "\t")
} else {
  message("QC note: neither ERBB2 nor GRB7 found; no QC summary written.")
}

# ===================== CONCORDANCE WITH TRANSCRIPT (optional if RNA DE exists) =====================
find_transcript_de <- function() {
  pat <- "DE_HER2pos_vs_HER2neg.*FDRlt0\\.05.*\\.(tsv|csv)$"
  cand <- unique(unlist(lapply(search_roots, function(root) {
    if (!dir.exists(root)) return(character())
    list.files(root, pattern = pat, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  })))
  if (!length(cand)) return(NA_character_)
  score <- function(p) {
    s <- 0
    s <- s + 3 * grepl("RNA|Transcript", p, ignore.case = TRUE)
    s <- s + 2 * grepl("Raj|RK|APOLLO4C", p, ignore.case = TRUE)
    s <- s + 1 * grepl("aim2_outputs", p, ignore.case = TRUE)
    s
  }
  cand[order(vapply(cand, score, numeric(1)), decreasing = TRUE)][1]
}

rna_de_path <- find_transcript_de()

if (is.character(rna_de_path) && file.exists(rna_de_path)) {
  message("Found transcript DE for concordance: ", rna_de_path)
  RNA <- tryCatch(fread(rna_de_path), error = function(e) NULL)
  
  if (!is.null(RNA) && nrow(RNA)) {
    cn <- names(RNA)
    id_col <- if ("feature_id" %in% cn) "feature_id" else if ("Gene" %in% cn) "Gene" else cn[1]
    lfcc <- cn[grepl("log2?FC", cn, ignore.case = TRUE)][1]
    if (is.na(lfcc)) lfcc <- cn[grepl("^logFC$", cn, ignore.case = TRUE)][1]
    
    setnames(RNA, id_col, "feature_id", skip_absent = TRUE)
    if (!is.na(lfcc)) setnames(RNA, lfcc, "log2FC_RNA", skip_absent = TRUE)
    
    RNA$feature_id <- toupper(trimws(as.character(RNA$feature_id)))
    RNA <- unique(RNA[, .(feature_id, log2FC_RNA)])
    
    PROT_core <- data.table(
      feature_id  = toupper(sig_FDR$feature_id),
      log2FC_PROT = sig_FDR$logFC
    )
    
    CC <- merge(PROT_core, RNA, by = "feature_id", all.x = TRUE, all.y = FALSE)
    CC[, concordant := !is.na(log2FC_RNA) & (log2FC_RNA * log2FC_PROT) > 0]
    
    fwrite(CC[!is.na(log2FC_RNA)], file.path(out_dir, "Concordant_PROT_vs_RNA.tsv"), sep = "\t")
    
    # FIXED: vector-safe concordance flag writing
    add_cc_flag <- function(df) {
      if (!nrow(df)) return(df)
      tmp <- merge(
        data.table(feature_id = toupper(df$feature_id)),
        CC[, .(feature_id, concordant)],
        by = "feature_id", all.x = TRUE
      )
      df$concordant_with_RNA <- ifelse(tmp$concordant == TRUE, "yes",
                                       ifelse(is.na(tmp$concordant), NA, "no"))
      df
    }
    
    if (file.exists(f_fdr)) {
      df <- fread(f_fdr, data.table = FALSE)
      df <- add_cc_flag(df)
      fwrite(df, f_fdr, sep = "\t")
    }
    if (file.exists(f_fc2)) {
      df <- fread(f_fc2, data.table = FALSE)
      df <- add_cc_flag(df)
      fwrite(df, f_fc2, sep = "\t")
    }
    
    CCp <- CC[!is.na(log2FC_RNA)]
    if (nrow(CCp) >= 3) {
      r <- suppressWarnings(cor(CCp$log2FC_RNA, CCp$log2FC_PROT, method = "pearson"))
      
      png(file.path(out_dir, "Concordant_scatter_log2FC_RNA_vs_PROT.png"), width = 1200, height = 1000, res = 150)
      par(mar = c(5,5,4.5,2))
      plot(CCp$log2FC_RNA, CCp$log2FC_PROT, pch = 21, bg = "white", col = "black",
           xlab = "log2FC (Transcript, HER2+ vs HER2−)",
           ylab = "log2FC (Proteome, HER2+ vs HER2−)",
           main = sprintf("Transcript–Proteome Concordance (n=%d, r=%.2f)", nrow(CCp), r))
      abline(lm(log2FC_PROT ~ log2FC_RNA, data = as.data.frame(CCp)), lty = 2, col = "grey40")
      grid(col = "grey85")
      dev.off()
    }
  }
} else {
  message("No transcript DE table found for concordance; skipped.")
}

# -------- VOLCANO (PNG-only; -log10(FDR)) --------
p05 <- -log10(0.05)

xv <- tt$logFC
yv <- -log10(pmax(tt$adj.P.Val, 1e-300))
labs <- tt$feature_id

is_sig_fdr <- tt$adj.P.Val < 0.05
is_fc1     <- abs(tt$logFC) >= 1

png(file.path(out_dir, "volcano.png"), width = 1600, height = 1100, res = 170)
par(mar = c(5,5,4.5,2))

plot(xv, yv,
     pch = 21,
     bg  = ifelse(is_sig_fdr, "black", "white"),
     col = "black",
     xlab = "log2FC (HER2pos - HER2neg)",
     ylab = expression(-log[10]("FDR")),
     main = cohort_id)

abline(v = c(-1, 1), lty = 2, col = "grey55")
abline(h = p05,      lty = 2, col = "grey55")
grid(col = "grey90")

force_genes <- c("ERBB2", "GRB7", "STARD3", "MIEN1")
idx_force <- which(labs %in% force_genes)
idx_core  <- which(is_sig_fdr & is_fc1)

max_core_labels <- 5
if (length(idx_core) > max_core_labels) {
  idx_core <- idx_core[order(tt$adj.P.Val[idx_core])][1:max_core_labels]
}

topN_extra <- 12
idx_top <- order(tt$adj.P.Val)
idx_top <- idx_top[!idx_top %in% c(idx_force, idx_core)]
idx_top <- head(idx_top, topN_extra)

idx_label <- unique(c(idx_force, idx_core, idx_top))

dx <- rep(c(-0.08, 0.08), length.out = length(idx_label))
dy <- rep(c(0.12, 0.18), length.out = length(idx_label))

text(xv[idx_label] + dx, yv[idx_label] + dy,
     labels = labs[idx_label],
     cex = 0.75,
     pos = ifelse(dx > 0, 4, 2),
     xpd = NA)

points(xv[idx_label], yv[idx_label],
       pch = 21,
       bg = ifelse(is_sig_fdr[idx_label] & is_fc1[idx_label], "black", "white"),
       col = "black")

dev.off()

# -------- RUN LOG --------
erbb2_line <- {
  i <- which(tt$feature_id == "ERBB2")[1]
  if (!is.na(i)) sprintf("ERBB2: log2FC=%.3f  FDR=%.3g", tt$logFC[i], tt$adj.P.Val[i]) else "ERBB2: NA"
}
grb7_line <- {
  i <- which(tt$feature_id == "GRB7")[1]
  if (!is.na(i)) sprintf("GRB7:  log2FC=%.3f  FDR=%.3g", tt$logFC[i], tt$adj.P.Val[i]) else "GRB7:  NA"
}

log_lines <- c(
  sprintf("RUN LOG — %s", cohort_id),
  sprintf("Proteome file: %s", basename(prot_file)),
  sprintf("Samples kept (after alignment/impute): %d (HER2pos=%d, HER2neg=%d)",
          ncol(M), sum(grp=="HER2pos"), sum(grp=="HER2neg")),
  sprintf("FDR<0.05: %d", nrow(sig_FDR)),
  sprintf("FDR<0.05 & |log2FC|>=1: %d", nrow(sig_FDR2)),
  "Meta:",
  paste0("  - ", meta_template_csv),
  paste0("  - ", meta_used_tsv),
  paste0("  - ", contrast_used_txt),
  "Outputs:",
  paste0("  - ", file.path(out_dir, "DE_full.tsv")),
  paste0("  - ", f_fdr),
  paste0("  - ", f_fc2),
  paste0("  - ", file.path(out_dir, "volcano.png")),
  paste0("  - ", file.path(out_dir, "qc_ERBB2_boxplot.png")),
  paste0("  - ", file.path(out_dir, "qc_GRB7_boxplot.png"))
)

if (is.character(rna_de_path) && file.exists(rna_de_path)) {
  log_lines <- c(
    log_lines,
    paste0("  - Concordance source (RNA): ", rna_de_path),
    paste0("  - ", file.path(out_dir, "Concordant_PROT_vs_RNA.tsv")),
    paste0("  - ", file.path(out_dir, "Concordant_scatter_log2FC_RNA_vs_PROT.png"))
  )
} else {
  log_lines <- c(log_lines, "  - Concordance: no transcript DE table found; skipped.")
}

log_lines <- c(log_lines, erbb2_line, grb7_line)
writeLines(log_lines, file.path(out_dir, "RUN_LOG.txt"))

message("DEP complete → ", out_dir)
# ================================= END ======================================
