# ============== run_Wolf_xCell.R ====================
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
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
labels_fp  <- file.path(base_dir, "labels_from_BP_subtype.csv")
expr_file  <- file.path(base_dir, "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt")

# >>> core HER2 activation features (Wolf RNA DE)
core_file  <- file.path(
  base_dir,
  "aim2_outputs",
  "GSE194040_ISPY2_mRNA_BPsubtype",
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

# ---- labels + HER2 groups ----------------------------------------
lab <- fread(labels_fp, data.table = FALSE)
id_col <- intersect(colnames(lab), c("SampleID","sample_id","Sample","ID","ResID","resid","ISPY2ResID"))
if (!length(id_col)) id_col <- colnames(lab)[1]
colnames(lab)[match(id_col[1], colnames(lab))] <- "SampleID_raw"

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
wanted    <- c("ImmuneScore", "StromaScore", "MicroenvironmentScore", "Fibroblasts", "Epithelial cells")
keep      <- intersect(wanted, sig_names)
if (!length(keep)) stop("No requested admixture signatures found in xCell results.")
message("Admixture signatures found: ", paste(keep, collapse = ", "))

subS <- S[keep, , drop = FALSE]

admix_dt <- data.table(
  SampleID_norm = colnames(subS),
  HER2          = as.character(grp)
)
for (sig in keep) admix_dt[[sig]] <- as.numeric(subS[sig, ])

admix_fp <- file.path(outdir, "Wolf_xCell_admixture_scores.tsv")
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
    title = "xCell admixture scores by HER2 group (Wolf)",
    x = "HER2 group",
    y = "xCell score"
  ) +
  theme(legend.position = "none")

ggsave(file.path(outdir, "Wolf_xCell_admixture_boxplots_by_HER2.pdf"),
       p_box, width = 8.2, height = 5.0)

# ---- Immune vs Stroma scatter (Wolf) ------------------------------
if (all(c("ImmuneScore","StromaScore") %in% colnames(admix_dt))) {
  p_scatter <- ggplot(admix_dt, aes(x = StromaScore, y = ImmuneScore, color = HER2)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_color_manual(values = COLS) +
    labs(
      title = "Immune vs Stroma admixture scores (Wolf)",
      x = "StromaScore",
      y = "ImmuneScore",
      color = "HER2"
    )
  ggsave(file.path(outdir, "Wolf_xCell_admixture_Immune_vs_Stroma_scatter.pdf"),
         p_scatter, width = 5.8, height = 4.6)
}

# ---- correlate core HER2 features with admixture ------------------
if (!file.exists(core_file)) {
  message("Core feature file not found: ", core_file, " — skipping correlations.")
} else if (!file.exists(expr_file)) {
  message("Expression file not found: ", expr_file, " — skipping correlations.")
} else {
  
  core <- fread(core_file)
  if (!"feature_id" %in% colnames(core)) stop("Core feature file must have a column named 'feature_id'.")
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
        message("No usable correlations computed — skipping correlation outputs.")
      } else {
        cor_dt <- rbindlist(cor_list)
        cor_dt[, FDR := p.adjust(p_value, method = "BH")]
        
        cor_fp <- file.path(outdir, "Wolf_xCell_admixture_correlations_core_features.tsv")
        fwrite(cor_dt, cor_fp, sep = "\t")
        message("Admixture–core feature correlations written to: ", cor_fp)
        
        # directional summary (strong positive vs strong negative)
        summary_dt <- cor_dt[, .(
          n_tests    = .N,
          n_pos      = sum(rho >  0.3 & FDR < 0.05),
          n_neg      = sum(rho < -0.3 & FDR < 0.05),
          n_strong   = sum(abs(rho) > 0.3 & FDR < 0.05),
          pct_pos    = 100 * sum(rho >  0.3 & FDR < 0.05) / .N,
          pct_neg    = 100 * sum(rho < -0.3 & FDR < 0.05) / .N,
          pct_strong = 100 * sum(abs(rho) > 0.3 & FDR < 0.05) / .N
        ), by = score_type]
        
        sum_fp <- file.path(outdir, "Wolf_xCell_admixture_correlations_summary_directional.tsv")
        fwrite(summary_dt, sum_fp, sep = "\t")
        message("Directional summary written to: ", sum_fp)
        
        # stacked barplot of % strong correlations by direction
        bar_dt <- rbind(
          summary_dt[, .(score_type, direction = "positive", pct = pct_pos)],
          summary_dt[, .(score_type, direction = "negative", pct = pct_neg)]
        )
        
        p_bar_dir <- ggplot(bar_dt, aes(x = score_type, y = pct, fill = direction)) +
          geom_col() +
          labs(
            title = "Core HER2 features linked to admixture (|rho| > 0.3, FDR < 0.05)",
            x = "xCell score",
            y = "% of core features",
            fill = "Direction"
          ) +
          theme(axis.text.x = element_text(angle = 30, hjust = 1))
        
        ggsave(file.path(outdir, "Wolf_xCell_admixture_correlations_barplot_directional.pdf"),
               p_bar_dir, width = 7.2, height = 4.6)
        
        # heatmap: top core genes vs xCell scores (rho)
        heat_dt <- cor_dt[, .(feature_id, score_type, rho, FDR)]
        top_genes <- heat_dt[, .(minFDR = suppressWarnings(min(FDR, na.rm = TRUE))), by = feature_id][order(minFDR)]
        top_genes <- top_genes[is.finite(minFDR)]
        top_genes <- head(top_genes$feature_id, 40)  # slide-safe size (adjust if needed)
        
        heat_sub <- heat_dt[feature_id %in% top_genes]
        mat_rho <- dcast(heat_sub, feature_id ~ score_type, value.var = "rho")
        rn <- mat_rho$feature_id
        mat <- as.matrix(mat_rho[, -1, drop = FALSE])
        rownames(mat) <- rn
        
        # cap extremes for readability
        mat_cap <- pmax(pmin(mat, 0.8), -0.8)
        
        pheatmap(
          mat_cap,
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          border_color = NA,
          main = "Core HER2 features vs xCell admixture (Spearman rho)",
          filename = file.path(outdir, "Wolf_xCell_core_features_vs_admixture_heatmap.pdf"),
          width = 7.5, height = 9
        )
        
        # example gene scatterplots (ERBB2/GRB7/STARD3/MIEN1 vs key scores)
        example_genes <- intersect(c("ERBB2","GRB7","STARD3","MIEN1"), rownames(E_sub))
        example_scores <- intersect(c("StromaScore","ImmuneScore","MicroenvironmentScore"), keep)
        
        for (g in example_genes) {
          df_g <- data.table(
            SampleID_norm = m$SampleID_norm,
            HER2 = admix_sub$HER2,
            expr = as.numeric(E_sub[g, ])
          )
          for (sc in example_scores) {
            df_g[[sc]] <- as.numeric(admix_sub[[sc]])
            p_ex <- ggplot(df_g, aes_string(x = sc, y = "expr", color = "HER2")) +
              geom_point(alpha = 0.6, size = 1.6) +
              scale_color_manual(values = COLS) +
              labs(
                title = paste0(g, " expression vs ", sc, " (Wolf)"),
                x = sc,
                y = paste0(g, " expression"),
                color = "HER2"
              )
            ggsave(file.path(outdir, paste0("Wolf_xCell_example_", g, "_vs_", sc, ".pdf")),
                   p_ex, width = 5.8, height = 4.4)
          }
        }
      }
    }
  }
}

message("Done. Outputs in: ", outdir)
# ==========================================================
