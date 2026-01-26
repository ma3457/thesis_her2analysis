# ================= run_Krug_proteome.R (PNG-only + Volcano = -log10(FDR)) =================
# Proteome DEP with LIMMA (HER2+ vs HER2−), dual HER2 definitions:
#   A) Clinical-only
#   B) Clinical + rescue (amplified/proteogenomic)
# Produces per tag:
#   - meta USED + contrast USED
#   - DE_full.tsv + filtered (FDR / FDR+FC2) + canonical copies
#   - volcano.png (y = -log10(adj.P.Val))
#   - MA.png
#   - QC_ERBB2/GRB7_boxplot.png
#   - RNA–PROT concordance (+ scatter.png)
#   - RUN_LOG.txt
# =========================================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(stringr)
})

# --------- INPUTS ----------
base_dir   <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/krug et al"
prot_file  <- file.path(base_dir, "HS_CPTAC_BRCA_2018_Proteome_Ratio_Norm_gene_Median.cct")
cli_file   <- file.path(base_dir, "HS_CPTAC_BRCA_2018_CLI.tsv")
out_prefix <- "Krug_CPTAC_BRCA_PROT"

stopifnot(file.exists(prot_file), file.exists(cli_file))

# --------- OUTPUT ROOT ----------
outputs_root <- file.path(base_dir, "aim2_outputs")
dir.create(outputs_root, showWarnings = FALSE, recursive = TRUE)
out_dir_tag <- function(tag) file.path(outputs_root, paste0(out_prefix, "_", tag))

# --------- SETTINGS ----------
min_per_group <- 3
make_concordance <- TRUE   # set FALSE if you want to skip RNA lookup

# ---- transcript DE lookup for concordance ----
search_roots <- c(base_dir, "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA")

find_transcript_de <- function(tag){
  pat <- "DE_HER2pos_vs_HER2neg.*FDRlt0\\.05.*\\.(tsv|csv)$"
  cand <- unique(unlist(lapply(search_roots, function(root){
    if (!dir.exists(root)) character()
    else list.files(root, pattern = pat, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  })))
  if (!length(cand)) return(NA_character_)
  
  score <- function(p){
    s <- 0
    s <- s + 3 * grepl("RNA|Transcript", p, ignore.case = TRUE)
    s <- s + 2 * grepl("Krug|CPTAC|BRCA", p, ignore.case = TRUE)
    s <- s + 1 * grepl("aim2_outputs", p, ignore.case = TRUE)
    s
  }
  cand[order(vapply(cand, score, numeric(1)), decreasing = TRUE)][1]
}

# --------- LOAD (Proteome: log2 TMT ratios) ----------
prot <- fread(prot_file, data.table = FALSE, check.names = FALSE)
rownames(prot) <- toupper(trimws(prot[[1]]))
prot[[1]] <- NULL
prot <- as.matrix(prot)
storage.mode(prot) <- "double"

cli <- fread(cli_file, data.table = FALSE, check.names = FALSE)
stopifnot("Sample.ID" %in% names(cli))

common <- intersect(colnames(prot), cli$Sample.ID)
prot <- prot[, common, drop = FALSE]
cli  <- cli[match(common, cli$Sample.ID), , drop = FALSE]
stopifnot(all(colnames(prot) == cli$Sample.ID))

# --------- HER2 GROUPS ----------
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

clinHER2 <- to_num(cli$Her2.Updated.Clinical.Status)  # 0/1/2/NA
ampHER2  <- if ("HER2.Amplified" %in% names(cli)) to_num(cli$HER2.Amplified) else rep(NA, nrow(cli))
protHER2 <- if ("ERBB2.Proteogenomic.Status" %in% names(cli)) to_num(cli$ERBB2.Proteogenomic.Status) else rep(NA, nrow(cli))

# A) Clinical-only: 1=pos, 0=neg, else NA
grpA <- ifelse(clinHER2 %in% c(0,1), ifelse(clinHER2==1,"HER2pos","HER2neg"), NA_character_)

# B) Clinical + rescue: inherit A; set missing to pos if amplified/proteogenomic pos
grpB <- grpA
rescue <- is.na(grpB) & (ampHER2 == 1 | protHER2 == 1)
grpB[rescue] <- "HER2pos"

# --------- QC PLOTS (PNG only) ----------
qc_box <- function(E, g_factor, tag, gene) {
  gene <- toupper(gene)
  if (!(gene %in% rownames(E))) return(invisible(FALSE))
  out_dir <- out_dir_tag(tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  fn_png <- file.path(out_dir, paste0("QC_", gene, "_boxplot.png"))
  png(fn_png, width = 1400, height = 1100, res = 160)
  par(mar = c(5,5,4.5,2))
  
  boxplot(split(as.numeric(E[gene, ]), g_factor),
          main = sprintf("%s QC — %s (n+=%d, n-=%d)",
                         gene, tag, sum(g_factor=="HER2pos"), sum(g_factor=="HER2neg")),
          xlab = "HER2 group", ylab = paste0(gene, " (log2 TMT)"),
          col = "grey90", border = "grey30")
  stripchart(split(as.numeric(E[gene, ]), g_factor), vertical = TRUE, method = "jitter",
             add = TRUE, pch = 21, bg = "black", col = "black")
  dev.off()
  TRUE
}

# --------- WRITER ----------
write_outputs <- function(res, tag, E_used, g_factor, full_fitAmean = NULL){
  
  out_dir <- out_dir_tag(tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # meta USED + contrast USED
  meta_used_tsv     <- file.path(out_dir, paste0(out_prefix, "_", tag, "_meta_USED.tsv"))
  contrast_used_txt <- file.path(out_dir, paste0(out_prefix, "_", tag, "_contrast_USED.txt"))
  
  meta_used <- data.table(
    SampleID = colnames(E_used),
    Group    = as.character(g_factor)
  )
  meta_used[, Her2_Updated_Clinical_Status := clinHER2[match(SampleID, cli$Sample.ID)]]
  meta_used[, HER2_Amplified               := ampHER2 [match(SampleID, cli$Sample.ID)]]
  meta_used[, ERBB2_Proteogenomic_Status   := protHER2[match(SampleID, cli$Sample.ID)]]
  fwrite(meta_used, meta_used_tsv, sep = "\t")
  
  writeLines(c(
    "Mode: HER2",
    "Contrast: HER2pos - HER2neg",
    paste0("Definition tag: ", tag),
    paste0("n_pos: ", sum(g_factor=="HER2pos")),
    paste0("n_neg: ", sum(g_factor=="HER2neg"))
  ), contrast_used_txt)
  
  # decorate results
  res$protein         <- toupper(res$protein)
  res$cohort_id       <- paste0(out_prefix, "_", tag)
  res$HER2_definition <- tag
  
  n_pos <- as.integer(sum(g_factor == "HER2pos"))
  n_neg <- as.integer(sum(g_factor == "HER2neg"))
  res$n_pos <- n_pos
  res$n_neg <- n_neg
  
  # Full DE
  full_path <- file.path(out_dir, "DE_full.tsv")
  fwrite(res, full_path, sep = "\t")
  
  # Filtered (FDR)
  sig <- res[is.finite(res$adj.P.Val) & res$adj.P.Val < 0.05,
             c("protein","log2FC_HER2pos_vs_HER2neg","AveExpr","t","P.Value","adj.P.Val","B",
               "cohort_id","HER2_definition","n_pos","n_neg")]
  
  f_fdr <- file.path(out_dir, paste0(out_prefix, "_", tag, "_DE_filtered_FDRlt0.05.tsv"))
  fwrite(sig, f_fdr, sep = "\t")
  
  sig_fc2 <- sig[is.finite(sig$log2FC_HER2pos_vs_HER2neg) & abs(sig$log2FC_HER2pos_vs_HER2neg) >= 1, ]
  f_fc2 <- file.path(out_dir, paste0(out_prefix, "_", tag, "_DE_filtered_FDRlt0.05_absLog2FCge1.tsv"))
  fwrite(sig_fc2, f_fc2, sep = "\t")
  
  # Volcano (PNG) — y = -log10(FDR)
  volc_path <- file.path(out_dir, "volcano.png")
  png(volc_path, width = 1400, height = 1100, res = 160)
  par(mar = c(5,5,4.5,2))
  
  # avoid Inf when FDR==0
  fdr_floor <- 1e-300
  y <- -log10(pmax(res$adj.P.Val, fdr_floor))
  x <- res$log2FC_HER2pos_vs_HER2neg
  
  plot(x, y,
       xlab = "log2FC (HER2pos - HER2neg)",
       ylab = "-log10(FDR)",
       main = paste0("Volcano — ", tag),
       pch = 21,
       bg  = ifelse(res$adj.P.Val < 0.05, "black", "white"),
       col = "black")
  
  abline(h = -log10(0.05), v = c(-1,1), lty = 2, col = "grey50")
  
  lab <- unique(c(head(res[order(res$adj.P.Val), "protein"], 10), "ERBB2","GRB7"))
  lab <- lab[lab %in% res$protein]
  sel <- res$protein %in% lab
  if (any(sel)) text(x[sel], y[sel], labels = res$protein[sel], pos = 4, cex = 0.8)
  dev.off()
  
  # MA (PNG)
  ma_path <- file.path(out_dir, "MA.png")
  png(ma_path, width = 1400, height = 1100, res = 160)
  par(mar = c(5,5,4.5,2))
  
  # Use AveExpr for MA x-axis (stable + already present in limma topTable)
  plot(res$AveExpr, res$log2FC_HER2pos_vs_HER2neg,
       xlab = "AveExpr", ylab = "log2FC",
       main = paste0("MA — ", tag),
       pch = 21, bg = "white", col = "black")
  abline(h = 0, lty = 2, col = "grey50")
  dev.off()
  
  # QC
  qc_box(E_used, g_factor, tag, "ERBB2")
  qc_box(E_used, g_factor, tag, "GRB7")
  
  # Concordance (optional)
  conc_tab  <- file.path(out_dir, paste0("Concordant_PROT_vs_RNA_", tag, ".tsv"))
  conc_plot <- file.path(out_dir, paste0("Concordant_scatter_log2FC_RNA_vs_PROT_", tag, ".png"))
  conc_log  <- file.path(out_dir, "CONCORDANCE_LOG.txt")
  
  if (isTRUE(make_concordance)) {
    rna_de_path <- find_transcript_de(tag)
    if (is.character(rna_de_path) && file.exists(rna_de_path)) {
      message("Found transcript DE for concordance: ", rna_de_path)
      RNA <- tryCatch(fread(rna_de_path), error = function(e) NULL)
      
      if (!is.null(RNA) && nrow(RNA)) {
        cn <- names(RNA)
        id_col <- if ("feature_id" %in% cn) "feature_id" else if ("Gene" %in% cn) "Gene" else cn[1]
        lfcc   <- cn[grepl("log2?FC|^logFC$", cn, ignore.case = TRUE)][1]
        
        setnames(RNA, id_col, "feature_id", skip_absent = TRUE)
        setnames(RNA, lfcc,   "log2FC_RNA",  skip_absent = TRUE)
        
        RNA[, feature_id := toupper(trimws(as.character(feature_id)))]
        RNA <- unique(RNA[, .(feature_id, log2FC_RNA)])
        
        PROT_core <- data.table(
          feature_id   = toupper(sig$protein),
          log2FC_PROT  = sig$log2FC_HER2pos_vs_HER2neg
        )
        
        CC <- merge(PROT_core, RNA, by = "feature_id", all.x = TRUE)
        CC[, concordant := !is.na(log2FC_RNA) & (log2FC_RNA * log2FC_PROT) > 0]
        fwrite(CC[!is.na(log2FC_RNA)], conc_tab, sep = "\t")
        
        CCp <- CC[!is.na(log2FC_RNA)]
        if (nrow(CCp) >= 3) {
          r <- suppressWarnings(cor(CCp$log2FC_RNA, CCp$log2FC_PROT, method = "pearson"))
          png(conc_plot, width = 1400, height = 1100, res = 160)
          par(mar = c(5,5,4.5,2))
          plot(CCp$log2FC_RNA, CCp$log2FC_PROT, pch = 21, bg = "white", col = "black",
               xlab = "log2FC (Transcript, HER2+ vs HER2−)",
               ylab = "log2FC (Proteome, HER2+ vs HER2−)",
               main = sprintf("Transcript–Proteome Concordance — %s (n=%d, r=%.2f)", tag, nrow(CCp), r))
          abline(lm(log2FC_PROT ~ log2FC_RNA, data = as.data.frame(CCp)), lty = 2, col = "grey40")
          grid(col = "grey85")
          dev.off()
        }
        
        writeLines(sprintf("Concordance source (RNA): %s", rna_de_path), conc_log)
      } else {
        writeLines("Transcript DE found but could not be read.", conc_log)
      }
    } else {
      writeLines("No transcript DE table found for concordance in configured search roots.", conc_log)
    }
  } else {
    writeLines("Concordance disabled (make_concordance=FALSE).", conc_log)
  }
  
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
    sprintf("Proteome file: %s", basename(prot_file)),
    sprintf("Samples kept: %d (HER2pos=%d, HER2neg=%d)", length(g_factor), n_pos, n_neg),
    sprintf("FDR<0.05: %d", nrow(sig)),
    sprintf("FDR<0.05 & |log2FC|>=1: %d", nrow(sig_fc2)),
    "Meta/contrast:",
    paste0("  - ", meta_used_tsv),
    paste0("  - ", contrast_used_txt),
    "Outputs:",
    paste0("  - ", full_path),
    paste0("  - ", f_fdr),
    paste0("  - ", f_fc2),
    paste0("  - ", volc_path),
    paste0("  - ", ma_path),
    paste0("  - ", file.path(out_dir, "QC_ERBB2_boxplot.png")),
    paste0("  - ", file.path(out_dir, "QC_GRB7_boxplot.png")),
    paste0("  - ", conc_log),
    erbb2_line,
    grb7_line
  )
  if (file.exists(conc_tab))  log_lines <- c(log_lines, paste0("  - ", conc_tab))
  if (file.exists(conc_plot)) log_lines <- c(log_lines, paste0("  - ", conc_plot))
  
  writeLines(log_lines, file.path(out_dir, "RUN_LOG.txt"))
  
  invisible(list(n_fdr = nrow(sig), n_fc2 = nrow(sig_fc2)))
}

# --------- LIMMA RUNNER ----------
run_dep <- function(E, group_vec, tag){
  
  keep <- !is.na(group_vec)
  E <- E[, keep, drop = FALSE]
  g <- factor(group_vec[keep], levels = c("HER2neg","HER2pos"))
  
  if (length(levels(g)) < 2 || min(table(g)) < min_per_group) {
    stop(sprintf("[%s] Need ≥%d per group. Counts: %s",
                 tag, min_per_group, paste(names(table(g)), table(g), collapse = ", ")))
  }
  
  # drop zero-variance proteins
  sds <- apply(E, 1, function(x) sd(x, na.rm = TRUE))
  sds[is.na(sds)] <- 0
  if (any(sds == 0)) E <- E[sds > 0, , drop = FALSE]
  
  design <- model.matrix(~ 0 + g)
  colnames(design) <- c("HER2neg","HER2pos")
  
  fit  <- lmFit(E, design)
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(HER2pos - HER2neg, levels = design)),
                 trend = TRUE, robust = TRUE)
  
  tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
  tt$protein <- rownames(tt)
  tt$log2FC_HER2pos_vs_HER2neg <- tt$logFC
  
  stats <- write_outputs(tt, tag, E, g)
  invisible(stats)
}

# --------- RUN ----------
statsA <- run_dep(prot, grpA, "A_ClinicalOnly")
statsB <- run_dep(prot, grpB, "B_ClinPlusRescue")

cat(sprintf("[DEP A_ClinicalOnly]    FDR<0.05=%d | FDR+FC2=%d\n", statsA$n_fdr, statsA$n_fc2))
cat(sprintf("[DEP B_ClinPlusRescue] FDR<0.05=%d | FDR+FC2=%d\n", statsB$n_fdr, statsB$n_fc2))
# ============================ END ============================
