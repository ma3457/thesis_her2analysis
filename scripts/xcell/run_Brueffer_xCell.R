# ============== run_Brueffer_xCell.R ===================
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

# ---- locate this script's folder -------------------------
script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) && file.exists(f)) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(".", winslash = "/")
})

base_dir  <- script_dir
res_dir   <- file.path(base_dir, "xcell_results")

# Brueffer HER2 mapping (from your prep folder)
labels_fp <- file.path(base_dir, "GSE81538_prep", "HER2_mapping_matched.csv")

# Brueffer RNA expression matrix (same base as DE script)
expr_file <- file.path(base_dir, "GSE81538_gene_expression_405_transformed.csv")

# >>> core HER2 activation features (Brueffer RNA DE; FC2 canonical from DE script)
core_file <- file.path(
  base_dir,
  "aim2_outputs",
  "GSE81538",
  "GSE81538_DE_filtered_FDRlt0.05_absLog2FCge1.tsv"
)

outdir <- file.path(base_dir, "aim2_outputs", "xcell_admixture")
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

# ---- xCell scores ------------------------------------
scores_fp <- list.files(res_dir, pattern = "_xCell_.*\\.txt$", full.names = TRUE)
scores_fp <- scores_fp[!grepl("RAW|pvals", scores_fp, ignore.case = TRUE)]
if (!length(scores_fp)) stop("No xCell score file found in xcell_results.")
if (length(scores_fp) > 1) stop("More than one xCell score file found; narrow pattern.")
message("Using xCell scores: ", basename(scores_fp))

dt <- fread(scores_fp, data.table = FALSE, check.names = FALSE)
rownames(dt) <- dt[[1]]
S <- as.matrix(dt[, -1, drop = FALSE])
colnames(S) <- norm_ids(colnames(S))

# ---- labels + HER2 groups (Brueffer mapping) ----------------------------------
lab <- fread(labels_fp, data.table = FALSE)

# Try to identify the sample column; fall back to first column if needed
id_col <- intersect(
  colnames(lab),
  c("sample_id","SampleID","SampleID_raw","Sample","sample","geo_accession","ID")
)
if (!length(id_col)) id_col <- colnames(lab)[1]
colnames(lab)[match(id_col[1], colnames(lab))] <- "SampleID_raw"

# ---- inspect label columns & pick HER2 status column --------------------------
message("Label file columns: ", paste(colnames(lab), collapse = ", "))

# 1) any column whose name contains 'her2' (case-insensitive)
cand_her2 <- grep("her2", colnames(lab), ignore.case = TRUE, value = TRUE)

# 2) if none, fall back to a generic 'group' column if it exists
if (!length(cand_her2) && "group" %in% colnames(lab)) {
  cand_her2 <- "group"
}

if (!length(cand_her2)) {
  stop(
    "Could not find a HER2 status column.\n",
    "Columns available: ", paste(colnames(lab), collapse = ", ")
  )
}

# use the first candidate and standardize the name
colnames(lab)[match(cand_her2[1], colnames(lab))] <- "her2_status"
message("Using HER2 status column: her2_status (from ", cand_her2[1], ")")

# Map to HER2pos / HER2neg (mirrors Wolf logic)
lab$HER2 <- ifelse(
  grepl("pos|3\\+|amp|amplif|overexpr|her2[-_ ]?high", lab$her2_status, ignore.case = TRUE),
  "HER2pos",
  ifelse(
    grepl("neg|0|1\\+|her2[-_ ]?low", lab$her2_status, ignore.case = TRUE),
    "HER2neg",
    NA_character_
  )
)

lab$SampleID_norm <- norm_ids(lab$SampleID_raw)

common <- intersect(colnames(S), lab$SampleID_norm)
if (!length(common)) stop("No overlap between xCell samples and labels.")
S   <- S[, common, drop = FALSE]
lab <- lab[match(common, lab$SampleID_norm), , drop = FALSE]
grp <- factor(lab$HER2, levels = c("HER2neg","HER2pos"))
stopifnot(all(!is.na(grp)))

# ---- admixture signatures -----------------------------------------------------
sig_names <- rownames(S)
wanted    <- c("ImmuneScore", "StromaScore", "MicroenvironmentScore", "Fibroblasts", "Epithelial cells")
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

admix_fp <- file.path(outdir, "Brueffer_xCell_admixture_scores.tsv")
fwrite(admix_dt, admix_fp, sep = "\t")
message("Admixture scores written to: ", admix_fp)

# ---- boxplots: admixture scores by HER2 ---------------------------------------
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
    title = "xCell admixture scores by HER2 group (Brueffer)",
    x = "HER2 group",
    y = "xCell score"
  ) +
  theme(legend.position = "none")

ggsave(file.path(outdir, "Brueffer_xCell_admixture_boxplots_by_HER2.pdf"),
       p_box, width = 8, height = 5)

# ---- scatter Immune vs Stroma ----------------------------------------
if (all(c("ImmuneScore","StromaScore") %in% colnames(admix_dt))) {
  p_scatter <- ggplot(admix_dt, aes(x = StromaScore, y = ImmuneScore, color = HER2)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_color_manual(values = COLS) +
    labs(
      title = "Immune vs Stroma admixture scores (Brueffer)",
      x = "StromaScore",
      y = "ImmuneScore",
      color = "HER2"
    )
  ggsave(file.path(outdir, "Brueffer_xCell_admixture_Immune_vs_Stroma_scatter.pdf"),
         p_scatter, width = 5.5, height = 4.5)
}

# ---- correlate core HER2 features with admixture -------------------
if (!file.exists(core_file)) {
  message("Core feature file not found: ", core_file,
          " — skipping correlations.")
} else if (!file.exists(expr_file)) {
  message("Expression matrix not found: ", expr_file,
          " — skipping correlations.")
} else {
  
  core <- fread(core_file)
  
  # Flexible feature ID detection
  id_col_core <- intersect(colnames(core), c("feature_id","gene_symbol","Gene","ID","Feature"))
  if (!length(id_col_core)) {
    stop(
      "Could not find a feature ID column in core_file.\n",
      "Columns available: ", paste(colnames(core), collapse = ", ")
    )
  }
  message("Using feature ID column in core_file: ", id_col_core[1])
  core_genes <- unique(core[[id_col_core[1]]])
  
  # RNA expression (same matrix as DE)
  x <- fread(file = expr_file, check.names = FALSE)
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
        
        cor_fp <- file.path(outdir, "Brueffer_xCell_admixture_correlations_core_features.tsv")
        fwrite(cor_dt, cor_fp, sep = "\t")
        message("Admixture–core feature correlations written to: ", cor_fp)
        
        # summary: % of features with |rho|>0.3 & FDR<0.05 per score
        summary_dt <- cor_dt[, .(
          n_tests    = .N,
          n_strong   = sum(abs(rho) > 0.3 & FDR < 0.05),
          pct_strong = 100 * sum(abs(rho) > 0.3 & FDR < 0.05) / .N
        ), by = score_type]
        
        sum_fp <- file.path(outdir, "Brueffer_xCell_admixture_correlations_summary.tsv")
        fwrite(summary_dt, sum_fp, sep = "\t")
        message("Admixture correlation summary written to: ", sum_fp)
        
        # barplot of % strong correlations per score
        p_bar <- ggplot(summary_dt, aes(x = score_type, y = pct_strong)) +
          geom_col(fill = "grey40") +
          ylim(0, max(summary_dt$pct_strong, 10)) +
          labs(
            title = "Core HER2 features linked to admixture scores (|rho| > 0.3 & FDR < 0.05)",
            x = "xCell admixture score",
            y = "% of core features"
          ) +
          theme(axis.text.x = element_text(angle = 30, hjust = 1))
        
        ggsave(file.path(outdir, "Brueffer_xCell_admixture_correlations_barplot.pdf"),
               p_bar, width = 6.5, height = 4.5)
      }
    }
  }
}

message("Done. Outputs in: ", outdir)
# =========================================================
