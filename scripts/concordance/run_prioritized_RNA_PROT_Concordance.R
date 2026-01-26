# ================= Prioritized (P76 / P25) RNA–Protein Concordance =================
# FINAL — robust column detection + explicit fread(file=) to avoid space-path issues
# Uses RNA consensus (median logFC across Wolf/Robinson/Brueffer), then tests protein support (RajKumar)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

# ---------------- CONFIG ----------------
FDR_MAX <- 0.05
FC_MIN  <- 1

OUT_BASE <- normalizePath(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_concordance_outputs/Prioritized_RNA_PROT_FINAL",
  mustWork = FALSE
)
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

# ---------------- INPUTS ----------------
P76_FILE <- normalizePath(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv",
  mustWork = FALSE
)
P25_FILE <- normalizePath(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_25_FDR0.05_LFC1_ALL3_with_stats_Wolf_Robinson_Brueffer.csv",
  mustWork = FALSE
)

WOLF_RNA <- normalizePath(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/GSE194040_ISPY2_mRNA_BPsubtype/ALL_HER2pos_vs_HER2neg_DE.tsv",
  mustWork = FALSE
)
ROBINSON_RNA <- normalizePath(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/robinson et al/aim2_outputs/GSE199633_Robinson2025/ALL_HER2pos_vs_HER2neg_DE.tsv",
  mustWork = FALSE
)
BRUEFFER_RNA <- normalizePath(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/brueffer et al/aim2_outputs/GSE81538/ALL_HER2pos_vs_HER2neg_DE.tsv",
  mustWork = FALSE
)

# UPDATED RajKumar PROT path (your DE_full.tsv)
PROT_FILE <- normalizePath(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/raj kumar et al/aim2_outputs/RajKumar_PROT_GLOBAL/DE_full.tsv",
  mustWork = FALSE
)

# ---------------- HELPERS ----------------
clean_gene <- function(x) toupper(str_trim(gsub("\\.\\d+$", "", as.character(x))))
sign_dir   <- function(x) ifelse(x > 0, "UP", ifelse(x < 0, "DOWN", NA))

canon <- function(x) {
  # lower + remove non-alphanum so adj.P.Val == adj_p_val == adjpval
  tolower(gsub("[^a-z0-9]+", "", x))
}

pick_col_ci <- function(nms, candidates) {
  nms_c  <- canon(nms)
  cand_c <- canon(candidates)
  hit <- match(cand_c, nms_c, nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit)) return(nms[hit[1]])
  NA_character_
}

stop_if_missing <- function(path, label) {
  if (is.na(path) || path == "" || !file.exists(path)) {
    stop("Missing file for ", label, ":\n", path)
  }
}

# ---------------- READ PRIORITIZED GENE SETS ----------------
read_gene_list <- function(path) {
  stop_if_missing(path, "prioritized set")
  dt <- fread(file = path)
  col <- names(dt)[1]
  unique(clean_gene(dt[[col]]))
}

# ---------------- READ RNA DE (robust) ----------------
read_rna_de <- function(path, label, mode) {
  stop_if_missing(path, paste0(label, " RNA DE"))
  dt <- fread(file = path)
  
  gene_col <- pick_col_ci(names(dt), c(
    "feature_id","gene_symbol","genesymbol","gene","genes","symbol","hgnc_symbol","Gene","SYMBOL","ID"
  ))
  fc_col <- pick_col_ci(names(dt), c(
    "logFC","log2FC","log2FoldChange","log2foldchange","LFC","log_fold_change","log2_ratio"
  ))
  fdr_col <- pick_col_ci(names(dt), c(
    "adj.P.Val","adjPVal","adjpval","adj_p_val","FDR","padj","qvalue","q_value","fdr"
  ))
  
  if (any(is.na(c(gene_col, fc_col, fdr_col)))) {
    cat("\n--- Column detection failed for:", label, "RNA ---\n")
    cat("File:", path, "\n")
    cat("Columns seen:\n")
    print(names(dt))
    stop("Cannot detect columns in ", label, " RNA DE.")
  }
  
  out <- dt[, .(
    gene  = clean_gene(get(gene_col)),
    logFC = as.numeric(get(fc_col)),
    fdr   = as.numeric(get(fdr_col))
  )]
  
  out <- out[!is.na(gene) & gene != "" & !is.na(logFC) & !is.na(fdr)]
  
  if (mode == "FDR_only") {
    out <- out[fdr <= FDR_MAX]
  } else {
    out <- out[fdr <= FDR_MAX & abs(logFC) >= FC_MIN]
  }
  
  out[, dir := sign_dir(logFC)]
  out <- out[!is.na(dir)]
  out[, cohort := label]
  
  unique(out)
}

# ---------------- READ PROTEIN DE (robust gene-level) ----------------
read_prot_de <- function(path) {
  stop_if_missing(path, "RajKumar PROT DE")
  dt <- fread(file = path)
  
  gene_col <- pick_col_ci(names(dt), c(
    "feature_id","gene_symbol","genesymbol","gene","genes","symbol","hgnc_symbol","Gene","SYMBOL","ID"
  ))
  fc_col <- pick_col_ci(names(dt), c(
    "logFC","log2FC","log2FoldChange","log2foldchange","LFC","log_fold_change","log2_ratio",
    "log2FC_HER2pos_vs_HER2neg","log2fc_her2pos_vs_her2neg"
  ))
  fdr_col <- pick_col_ci(names(dt), c(
    "adj.P.Val","adjPVal","adjpval","adj_p_val","FDR","padj","qvalue","q_value","fdr"
  ))
  
  if (any(is.na(c(gene_col, fc_col, fdr_col)))) {
    cat("\n--- Column detection failed for: RajKumar PROT ---\n")
    cat("File:", path, "\n")
    cat("Columns seen:\n")
    print(names(dt))
    stop("Cannot detect columns in RajKumar PROT DE.")
  }
  
  out <- dt[, .(
    gene = clean_gene(get(gene_col)),
    prot_logFC = as.numeric(get(fc_col)),
    prot_fdr   = as.numeric(get(fdr_col))
  )]
  
  out <- out[!is.na(gene) & gene != "" & !is.na(prot_logFC) & !is.na(prot_fdr)]
  out <- out[prot_fdr <= FDR_MAX]
  out[, prot_dir := sign_dir(prot_logFC)]
  out <- out[!is.na(prot_dir)]
  unique(out)
}

# ---------------- MAIN RUNNER ----------------
analyze_set <- function(set_name, genes, rna_all, rna_cons, prot, outdir, mode) {
  
  rna_set  <- rna_cons[gene %in% genes]
  prot_set <- prot[gene %in% genes]
  
  # consensus vs protein overlap
  m <- merge(rna_set, prot_set, by = "gene")
  m[, concordant := (rna_dir == prot_dir)]
  
  # save tables
  fwrite(m[order(-abs(prot_logFC))],
         file.path(outdir, paste0(set_name, "_overlap_table.tsv")),
         sep = "\t")
  
  summary <- data.table(
    set = set_name,
    mode = mode,
    genes_in_set = length(genes),
    genes_in_rna_consensus = nrow(rna_set),
    genes_in_proteomics_sig = nrow(prot_set),
    overlap_genes = nrow(m),
    concordant_genes = if (nrow(m)) sum(m$concordant) else 0,
    concordant_pct = if (nrow(m)) round(100 * sum(m$concordant) / nrow(m), 1) else NA_real_
  )
  fwrite(summary, file.path(outdir, paste0(set_name, "_summary.tsv")), sep = "\t")
  
  # OPTIONAL: show how many RNA cohorts support each overlapped gene (at RNA level)
  if (nrow(m) > 0) {
    rna_support <- rna_all[gene %in% m$gene, .(rna_support_n = uniqueN(cohort)), by = gene]
    m2 <- merge(m, rna_support, by = "gene", all.x = TRUE)
    fwrite(m2[order(-rna_support_n, -abs(prot_logFC))],
           file.path(outdir, paste0(set_name, "_overlap_with_RNA_support.tsv")),
           sep = "\t")
  }
  
  # plot
  if (nrow(m) > 0) {
    label_genes <- intersect(c("ERBB2","GRB7","STARD3","PGAP3","MIEN1","CDK12","CCNK","ERBB4"), m$gene)
    
    p <- ggplot(m, aes(rna_logFC, prot_logFC)) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
      geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
      geom_point(aes(shape = concordant), size = 3) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
      theme_bw(base_size = 14) +
      labs(
        title = paste0(set_name, " RNA–Protein Concordance"),
        subtitle = paste0(
          "Mode: ", mode,
          " | overlap=", nrow(m),
          " | concordant=", sum(m$concordant)
        ),
        x = "Consensus RNA log2FC (median across Wolf/Robinson/Brueffer)",
        y = "RajKumar protein log2FC (HER2+ vs HER2−)"
      )
    
    if (length(label_genes)) {
      p <- p + geom_text(
        data = m[gene %in% label_genes],
        aes(label = gene),
        nudge_x = 0.05, nudge_y = 0.05,
        size = 3.5
      )
    }
    
    ggsave(file.path(outdir, paste0(set_name, "_scatter.png")),
           p, width = 8.5, height = 5.5, dpi = 200)
    ggsave(file.path(outdir, paste0(set_name, "_scatter.pdf")),
           p, width = 8.5, height = 5.5)
  }
  
  invisible(summary)
}

run_mode <- function(mode) {
  
  outdir <- file.path(OUT_BASE, mode)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # RNA (3 cohorts)
  rna_all <- rbindlist(list(
    read_rna_de(WOLF_RNA, "Wolf", mode),
    read_rna_de(ROBINSON_RNA, "Robinson", mode),
    read_rna_de(BRUEFFER_RNA, "Brueffer", mode)
  ), fill = TRUE)
  
  fwrite(rna_all, file.path(outdir, "RNA_significant_all_cohorts.tsv"), sep = "\t")
  
  # RNA consensus (median across cohorts)
  rna_cons <- rna_all[, .(
    rna_logFC = median(logFC, na.rm = TRUE),
    rna_fdr_min = min(fdr, na.rm = TRUE),
    n_cohorts = uniqueN(cohort)
  ), by = gene]
  rna_cons[, rna_dir := sign_dir(rna_logFC)]
  rna_cons <- rna_cons[!is.na(rna_dir)]
  fwrite(rna_cons[order(-abs(rna_logFC))],
         file.path(outdir, "RNA_consensus_median.tsv"),
         sep = "\t")
  
  # protein
  prot <- read_prot_de(PROT_FILE)
  fwrite(prot[order(-abs(prot_logFC))],
         file.path(outdir, "PROT_significant.tsv"),
         sep = "\t")
  
  # gene sets
  P76 <- read_gene_list(P76_FILE)
  P25 <- read_gene_list(P25_FILE)
  
  # analyze
  s1 <- analyze_set("P76", P76, rna_all, rna_cons, prot, outdir, mode)
  s2 <- analyze_set("P25", P25, rna_all, rna_cons, prot, outdir, mode)
  
  # combined summary
  fwrite(rbindlist(list(s1, s2), fill = TRUE),
         file.path(outdir, "SUMMARY_P76_P25.tsv"),
         sep = "\t")
  
  cat("\nDONE:", mode, "->", outdir, "\n")
}

# ---------------- RUN BOTH MODES ----------------
run_mode("FDR_only")
run_mode("FDR_FC")

cat("\nALL DONE -> ", OUT_BASE, "\n", sep = "")
