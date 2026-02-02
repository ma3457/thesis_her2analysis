# ============== run_Robinson_xCell.R ===================
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

# locate this script's folder
script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) && file.exists(f)) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(".", winslash = "/")
})

base_dir   <- script_dir
res_dir    <- file.path(base_dir, "xcell_results")

# ---- metadata + expression + core features for ROBINSON ----------
labels_fp  <- file.path(
  base_dir,
  "GSE199633_Robinson2025_prep",
  "HER2_mapping_matched.csv"
)

expr_file  <- file.path(base_dir, "GSE199633_expr_for_pipeline.tsv")

# >>> core HER2 activation features (Robinson RNA DE)
core_file  <- file.path(
  base_dir,
  "aim2_outputs",
  "GSE199633_Robinson2025",
  "FILT_HER2pos_vs_HER2neg_DE_FDR0.05_LFC1.tsv"
)

outdir     <- file.path(base_dir, "aim2_outputs", "xcell_admixture")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

COLS <- c(HER2neg = "#1f77b4", HER2pos = "#d62728")
theme_set(theme_bw(base_size = 11))

norm_ids <- function(v) {
  v <- as.character(v)
  v <- sub("^X", "", v)
  v <- sub("\\.GPL\\d+$", "", v)
  v <- gsub("[^A-Za-z0-9]", "", v)
  toupper(v)
}

# ---- xCell scores -------------------------------------------------
scores_fp <- list.files(res_dir, pattern = "_xCell_.*\\.txt$", full.names = TRUE)
scores_fp <- scores_fp[!grepl("RAW|pvals", scores_fp, ignore.case = TRUE)]
if (!length(scores_fp)) stop("No xCell score file found in xcell_results.")
if (length(scores_fp) > 1) stop("More than one xCell score file found; narrow pattern.")
message("Using xCell scores: ", basename(scores_fp))

dt <- fread(scores_fp, data.table = FALSE, check.names = FALSE)
rownames(dt) <- dt[[1]]
S <- as.matrix(dt[, -1, drop = FALSE])
colnames(S) <- norm_ids(colnames(S))

# ---- labels + HER2 groups (ROBINSON) ------------------------------
lab <- fread(labels_fp, data.table = FALSE)
cn  <- tolower(colnames(lab))
colnames(lab) <- cn

id_candidates  <- c("sampleid","sample_id","sample","id","resid",
                    "geo_accession","gsm")
grp_candidates <- c("group","her2_status","her2","her2_simple","her2_status_simple")

id_col  <- intersect(id_candidates,  cn)[1]
grp_col <- intersect(grp_candidates, cn)[1]

if (is.na(id_col) || is.na(grp_col)) {
  stop("Could not find sample/group columns in labels file. Have columns: ",
       paste(colnames(lab), collapse = ", "))
}

lab$SampleID_raw <- lab[[id_col]]
lab$group        <- lab[[grp_col]]

stopifnot("group" %in% colnames(lab))

lab$HER2 <- ifelse(grepl("pos|\\+$", lab$group, ignore.case = TRUE), "HER2pos",
                   ifelse(grepl("neg|-$", lab$group, ignore.case = TRUE), "HER2neg", NA_character_))
lab$SampleID_norm <- norm_ids(lab$SampleID_raw)

common <- intersect(colnames(S), lab$SampleID_norm)
if (!length(common)) stop("No overlap between xCell samples and labels.")
S   <- S[, common, drop = FALSE]
lab <- lab[match(common, lab$SampleID_norm), , drop = FALSE]
grp <- factor(lab$HER2, levels = c("HER2neg","HER2pos"))
stopifnot(all(!is.na(grp)))

# ---- admixture signatures ----------------------------------------
sig_names <- rownames(S)
wanted    <- c("ImmuneScore", "StromaScore", "MicroenvironmentScore",
               "Fibroblasts", "Epithelial cells")
keep      <- intersect(wanted, sig_names)
if (!length(keep)) stop("No requested admixture signatures found in xCell results.")
message("Admixture signatures found: ", paste(keep, collapse = ", "))

subS <- S[keep, , drop = FALSE]

admix_dt <- data.table(
  SampleID_norm = colnames(subS),
  HER2          = as.character(grp)
)
for (sig in keep) {
  admix_dt[[sig]] <- as.numeric(subS[sig, ])
}

admix_fp <- file.path(outdir, "Robinson_xCell_admixture_scores.tsv")
fwrite(admix_dt, admix_fp, sep = "\t")
message("Admixture scores written to: ", admix_fp)

# ---- boxplots: admixture scores by HER2 ---------------------------
long_dt <- melt(
  admix_dt,
  id.vars = c("SampleID_norm", "HER2"),
  variable.name = "score_type",
  value.name = "score"
)

p_box <- ggplot(long_dt, aes(x = HER2, y = score, fill = HER2)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.9) +
  scale_fill_manual(values = COLS) +
  facet_wrap(~ score_type, scales = "free_y") +
  labs(
    title = "xCell admixture scores by HER2 group (Robinson)",
    x = "HER2 group",
    y = "xCell score"
  ) +
  theme(legend.position = "none")

ggsave(file.path(outdir, "Robinson_xCell_admixture_boxplots_by_HER2.pdf"),
       p_box, width = 8, height = 5)

# optional scatter Immune vs Stroma
if (all(c("ImmuneScore","StromaScore") %in% colnames(admix_dt))) {
  p_scatter <- ggplot(admix_dt, aes(x = StromaScore, y = ImmuneScore, color = HER2)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_color_manual(values = COLS) +
    labs(
      title = "Immune vs Stroma admixture scores (Robinson)",
      x = "StromaScore",
      y = "ImmuneScore",
      color = "HER2"
    )
  ggsave(file.path(outdir, "Robinson_xCell_admixture_Immune_vs_Stroma_scatter.pdf"),
         p_scatter, width = 5.5, height = 4.5)
}

# ---- correlate core HER2 features with admixture ------
if (!file.exists(core_file)) {
  message("Core feature file not found: ", core_file,
          " — skipping correlations.")
} else {
  
  core <- fread(core_file)
  if (!"feature_id" %in% colnames(core)) {
    stop("Core feature file must have a column named 'feature_id'.")
  }
  core_genes <- unique(core$feature_id)
  
  # RNA expression
  x <- fread(expr_file, check.names = FALSE)
  if (is.na(names(x)[1]) || names(x)[1] == "") setnames(x, 1, "gene_symbol")
  if (!names(x)[1] %in% c("gene_symbol","Gene","feature_id","ID","Feature"))
    setnames(x, 1, "gene_symbol")
  x <- x[!duplicated(x$gene_symbol), ]
  
  samp_cols_expr <- setdiff(names(x), "gene_symbol")
  for (nm in samp_cols_expr) {
    if (!is.numeric(x[[nm]])) {
      v <- as.character(x[[nm]]); v[is.na(v)] <- ""
      suppressWarnings(num <- as.numeric(v))
      if (sum(nchar(v) > 0 & is.na(num)) == 0) x[[nm]] <- num
    }
  }
  E <- as.matrix(x[, ..samp_cols_expr])
  rownames(E) <- x$gene_symbol
  storage.mode(E) <- "numeric"
  
  expr_ids_norm <- norm_ids(colnames(E))
  expr_map <- data.table(expr_col = colnames(E), SampleID_norm = expr_ids_norm)
  
  # align expression to admixture samples
  m <- merge(admix_dt[, .(SampleID_norm)], expr_map, by = "SampleID_norm")
  if (!nrow(m)) {
    message("No overlap between expression and xCell admixture samples — skipping correlations.")
  } else {
    
    E_sub <- E[, m$expr_col, drop = FALSE]
    stopifnot(ncol(E_sub) == nrow(m))
    admix_sub <- admix_dt[match(m$SampleID_norm, admix_dt$SampleID_norm), ]
    
    genes_use <- intersect(core_genes, rownames(E_sub))
    if (!length(genes_use)) {
      message("No overlap between core features and RNA matrix — skipping correlations.")
    } else {
      message("Core genes with RNA data: ", length(genes_use))
      
      score_cols <- keep
      cor_list <- list()
      
      for (g in genes_use) {
        expr_vec <- as.numeric(E_sub[g, ])
        for (sc in score_cols) {
          score_vec <- as.numeric(admix_sub[[sc]])
          if (all(is.na(score_vec))) next
          ok <- is.finite(expr_vec) & is.finite(score_vec)
          if (sum(ok) < 10) next
          ct <- suppressWarnings(cor.test(expr_vec[ok], score_vec[ok], method = "spearman"))
          cor_list[[length(cor_list) + 1L]] <- data.table(
            feature_id = g,
            score_type = sc,
            rho        = unname(ct$estimate),
            p_value    = ct$p.value,
            n_used     = sum(ok)
          )
        }
      }
      
      if (!length(cor_list)) {
        message("No usable correlations computed — skipping correlation plots.")
      } else {
        cor_dt <- rbindlist(cor_list)
        cor_dt[, FDR := p.adjust(p_value, method = "BH")]
        
        cor_fp <- file.path(outdir, "Robinson_xCell_admixture_correlations_core_features.tsv")
        fwrite(cor_dt, cor_fp, sep = "\t")
        message("Admixture–core feature correlations written to: ", cor_fp)
        
        # summary: % of features with |rho|>0.3 & FDR<0.05 per score
        summary_dt <- cor_dt[, .(
          n_tests      = .N,
          n_strong     = sum(abs(rho) > 0.3 & FDR < 0.05),
          pct_strong   = 100 * sum(abs(rho) > 0.3 & FDR < 0.05) / .N
        ), by = score_type]
        
        sum_fp <- file.path(outdir, "Robinson_xCell_admixture_correlations_summary.tsv")
        fwrite(summary_dt, sum_fp, sep = "\t")
        message("Admixture correlation summary written to: ", sum_fp)
        
        # barplot of % strong correlations per score
        p_bar <- ggplot(summary_dt, aes(x = score_type, y = pct_strong)) +
          geom_col(fill = "grey40") +
          ylim(0, max(summary_dt$pct_strong, 10)) +
          labs(
            title = "Core HER2 features linked to admixture scores (Robinson, |rho| > 0.3 & FDR < 0.05)",
            x = "xCell admixture score",
            y = "% of core features"
          ) +
          theme(axis.text.x = element_text(angle = 30, hjust = 1))
        
        ggsave(file.path(outdir, "Robinson_xCell_admixture_correlations_barplot.pdf"),
               p_bar, width = 6.5, height = 4.5)
      }
    }
  }
}

message("Done. Outputs in: ", outdir)
# ==========================================================
