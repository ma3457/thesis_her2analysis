# ================= run_prioritized_xCell_admixture.R  =================
# Phase 4: quantify admixture effects on PRIORITIZED HER2 feature sets (P76, P25)
# Cohorts: Wolf, Robinson, Brueffer
# Outputs: correlation tables (ALL, FDR-only, FDR+|rho|), directional summaries, slide-ready barplots,
#          plus immune vs stroma scatter per cohort (xCell)
# Added: cross-cohort SUMMARY table + one combined summary barplot

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# --------------------------- CONFIG ---------------------------------

# Prioritized sets (provided)
P76_FILE <- normalizePath(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv",
  mustWork = TRUE
)
P25_FILE <- normalizePath(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_25_FDR0.05_LFC1_ALL3_with_stats_Wolf_Robinson_Brueffer.csv",
  mustWork = TRUE
)

# Output root
OUT_DIR <- normalizePath(
  "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_outputs/xcell_admixture_prioritized",
  mustWork = FALSE
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# thresholds
FDR_MAX <- 0.05
RHO_MIN <- 0.30

# xCell scores to test (will auto-intersect with what's available)
WANTED_SCORES <- c("ImmuneScore", "StromaScore", "MicroenvironmentScore", "Fibroblasts", "Epithelial cells")

# colors
COLS <- c(HER2neg = "#1f77b4", HER2pos = "#d62728")

# Cohort-specific inputs
COHORTS <- list(
  Wolf = list(
    cohort_name = "Wolf",
    # xCell score file (the non-RAW, non-pvals one)
    xcell_scores_fp = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/xcell_results", # folder; script will auto-pick file
    # labels file used previously (has group column)
    labels_fp = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/labels_from_BP_subtype.csv",
    # expression matrix
    expr_fp   = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt",
    labels_type = "wolf"  # custom parsing
  ),
  Robinson = list(
    cohort_name = "Robinson",
    xcell_scores_fp = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/robinson et al/xcell_results/xCell_GSE199633_expr_for_pipeline_xCell_1812102625.txt",
    labels_fp = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/robinson et al/GSE199633_Robinson2025_prep/HER2_mapping_matched.csv",
    expr_fp   = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/robinson et al/GSE199633_expr_for_pipeline.tsv",
    labels_type = "matched"  # sample + HER2_status
  ),
  Brueffer = list(
    cohort_name = "Brueffer",
    xcell_scores_fp = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/brueffer et al/xcell_results/xCell_GSE81538_gene_expression_405_transformed_xCell_1848102625.txt",
    labels_fp = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/brueffer et al/GSE81538_prep/HER2_mapping_matched.csv",
    expr_fp   = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/brueffer et al/GSE81538_gene_expression_405_transformed.csv",
    labels_type = "matched"
  )
)

# --------------------------- HELPERS --------------------------------

clean_gene <- function(x) toupper(str_trim(gsub("\\.\\d+$", "", as.character(x))))

norm_ids <- function(v) {
  v <- as.character(v)
  v <- sub("^X", "", v)
  v <- sub("\\.GPL\\d+$", "", v)
  v <- gsub("[^A-Za-z0-9]", "", v)
  toupper(v)
}

# read prioritized gene lists
get_prioritized_genes <- function(csv_fp) {
  dt <- fread(file = csv_fp)
  gene_col <- intersect(names(dt), c("gene","gene_symbol","Gene","SYMBOL","hgnc_symbol"))
  if (!length(gene_col)) stop("Cannot find gene column in: ", csv_fp)
  unique(clean_gene(dt[[gene_col[1]]]))
}

# read xCell score matrix: rows=scores, cols=samples
read_xcell_scores <- function(path_or_file) {
  # allow folder (Wolf) or file (Robinson/Brueffer)
  if (dir.exists(path_or_file)) {
    fp <- list.files(path_or_file, pattern = "_xCell_.*\\.txt$", full.names = TRUE)
    fp <- fp[!grepl("RAW|pvals", fp, ignore.case = TRUE)]
    if (!length(fp)) stop("No xCell score file found in: ", path_or_file)
    if (length(fp) > 1) stop("Multiple xCell score files found in: ", path_or_file, " (narrow pattern).")
    scores_fp <- fp[1]
  } else {
    scores_fp <- path_or_file
  }
  
  dt <- fread(file = scores_fp, data.table = FALSE, check.names = FALSE)
  rownames(dt) <- dt[[1]]
  S <- as.matrix(dt[, -1, drop = FALSE])
  colnames(S) <- norm_ids(colnames(S))
  list(S = S, scores_fp = scores_fp)
}

# read labels
read_labels <- function(labels_fp, labels_type) {
  lab <- fread(file = labels_fp, data.table = FALSE)
  
  if (labels_type == "wolf") {
    # Wolf file uses a "group" column with pos/neg embedded
    if (!"group" %in% colnames(lab)) stop("Wolf labels must contain column 'group': ", labels_fp)
    
    id_col <- intersect(colnames(lab), c("SampleID","sample_id","Sample","ID","ResID","resid","ISPY2ResID"))
    if (!length(id_col)) id_col <- colnames(lab)[1]
    colnames(lab)[match(id_col[1], colnames(lab))] <- "SampleID_raw"
    
    g <- as.character(lab$group)
    her2 <- ifelse(grepl("pos|\\+$", g, ignore.case = TRUE), "HER2pos",
                   ifelse(grepl("neg|-$", g, ignore.case = TRUE), "HER2neg", NA_character_))
    lab$HER2 <- her2
    lab$SampleID_norm <- norm_ids(lab$SampleID_raw)
    lab <- lab[!is.na(lab$HER2), , drop = FALSE]
    return(lab[, c("SampleID_norm","HER2"), drop = FALSE])
  }
  
  if (labels_type == "matched") {
    # matched files have: sample, HER2_status
    req <- c("sample","HER2_status")
    if (!all(req %in% colnames(lab))) stop("Matched labels must have columns: sample, HER2_status: ", labels_fp)
    
    lab$SampleID_norm <- norm_ids(lab$sample)
    hs <- as.character(lab$HER2_status)
    her2 <- ifelse(grepl("^pos|positive$", hs, ignore.case = TRUE), "HER2pos",
                   ifelse(grepl("^neg|negative$", hs, ignore.case = TRUE), "HER2neg", NA_character_))
    lab$HER2 <- her2
    lab <- lab[!is.na(lab$HER2), , drop = FALSE]
    return(lab[, c("SampleID_norm","HER2"), drop = FALSE])
  }
  
  stop("Unknown labels_type: ", labels_type)
}

# expression matrix reader: first col=gene, rest=sample columns
read_expression_gene_matrix <- function(expr_fp) {
  x <- fread(file = expr_fp, check.names = FALSE)
  if (is.na(names(x)[1]) || names(x)[1] == "") setnames(x, 1, "gene_symbol")
  setnames(x, 1, "gene_symbol")
  
  x$gene_symbol <- clean_gene(x$gene_symbol)
  x <- x[!duplicated(x$gene_symbol), ]
  
  samp_cols <- setdiff(names(x), "gene_symbol")
  
  # coerce numeric where safe
  for (nm in samp_cols) {
    if (!is.numeric(x[[nm]])) {
      v <- as.character(x[[nm]]); v[is.na(v)] <- ""
      suppressWarnings(num <- as.numeric(v))
      if (sum(nchar(v) > 0 & is.na(num)) == 0) x[[nm]] <- num
    }
  }
  
  E <- as.matrix(x[, ..samp_cols])
  rownames(E) <- x$gene_symbol
  storage.mode(E) <- "numeric"
  
  expr_map <- data.table(expr_col = colnames(E), SampleID_norm = norm_ids(colnames(E)))
  list(E = E, expr_map = expr_map)
}

# ---------------------- LOAD PRIORITIZED SETS -----------------------

p76_genes <- get_prioritized_genes(P76_FILE)
p25_genes <- get_prioritized_genes(P25_FILE)

cat("P76 genes:", length(p76_genes), "\n")
cat("P25 genes:", length(p25_genes), "\n")

# ---------------------- MAIN PER COHORT -----------------------------

all_cor <- list()

for (nm in names(COHORTS)) {
  
  cfg <- COHORTS[[nm]]
  cohort <- cfg$cohort_name
  cat("\n--- Running cohort:", cohort, "---\n")
  
  cohort_out <- file.path(OUT_DIR, cohort)
  dir.create(cohort_out, recursive = TRUE, showWarnings = FALSE)
  
  # xCell
  sx <- read_xcell_scores(cfg$xcell_scores_fp)
  S <- sx$S
  cat("xCell scores file: ", sx$scores_fp, "\n", sep="")
  
  # labels
  lab <- read_labels(cfg$labels_fp, cfg$labels_type)
  
  # align xCell samples to labels
  common <- intersect(colnames(S), lab$SampleID_norm)
  if (!length(common)) stop("No overlap between xCell samples and labels for cohort: ", cohort)
  
  S <- S[, common, drop = FALSE]
  lab <- lab[match(common, lab$SampleID_norm), , drop = FALSE]
  grp <- factor(lab$HER2, levels = c("HER2neg","HER2pos"))
  
  # select scores available
  keep <- intersect(WANTED_SCORES, rownames(S))
  if (!length(keep)) stop("None of WANTED_SCORES found in xCell output for: ", cohort)
  
  subS <- S[keep, , drop = FALSE]
  admix_dt <- data.table(
    SampleID_norm = colnames(subS),
    HER2          = as.character(grp)
  )
  for (sig in keep) admix_dt[[sig]] <- as.numeric(subS[sig, ])
  
  fwrite(admix_dt, file.path(cohort_out, paste0(cohort, "_xCell_admixture_scores.tsv")), sep = "\t")
  
  # scatter (for your 3-panel slide)
  if (all(c("ImmuneScore","StromaScore") %in% colnames(admix_dt))) {
    p_sc <- ggplot(admix_dt, aes(x = StromaScore, y = ImmuneScore, color = HER2)) +
      geom_point(alpha = 0.6, size = 1.6) +
      scale_color_manual(values = COLS) +
      labs(
        title = paste0("Immune vs Stroma admixture scores (", cohort, ")"),
        x = "StromaScore",
        y = "ImmuneScore",
        color = "HER2"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(cohort_out, paste0(cohort, "_Immune_vs_Stroma_scatter.pdf")),
           p_sc, width = 5.8, height = 4.6)
  }
  
  # expression
  ex <- read_expression_gene_matrix(cfg$expr_fp)
  E <- ex$E
  expr_map <- ex$expr_map
  
  m <- merge(admix_dt[, .(SampleID_norm)], expr_map, by = "SampleID_norm")
  if (!nrow(m)) stop("No overlap between expression samples and xCell samples for: ", cohort)
  
  E_sub <- E[, m$expr_col, drop = FALSE]
  admix_sub <- admix_dt[match(m$SampleID_norm, admix_dt$SampleID_norm), ]
  
  correlate_set <- function(genes, set_name) {
    genes_use <- intersect(genes, rownames(E_sub))
    if (!length(genes_use)) {
      warning("No overlap between ", set_name, " and expression genes for cohort: ", cohort)
      return(NULL)
    }
    
    cor_list <- list()
    for (g in genes_use) {
      expr_vec <- as.numeric(E_sub[g, ])
      for (sc in keep) {
        score_vec <- as.numeric(admix_sub[[sc]])
        ok <- is.finite(expr_vec) & is.finite(score_vec)
        if (sum(ok) < 10) next
        ct <- suppressWarnings(cor.test(expr_vec[ok], score_vec[ok], method = "spearman"))
        cor_list[[length(cor_list) + 1L]] <- data.table(
          cohort     = cohort,
          set        = set_name,
          gene       = g,
          score_type = sc,
          rho        = unname(ct$estimate),
          p_value    = ct$p.value,
          n_used     = sum(ok)
        )
      }
    }
    if (!length(cor_list)) return(NULL)
    
    cor_dt <- rbindlist(cor_list)
    # adjust within cohort+set for each score_type (keeps logic consistent across cohorts)
    cor_dt[, FDR := p.adjust(p_value, method = "BH"), by = .(score_type)]
    
    # output tables
    fwrite(cor_dt, file.path(cohort_out, paste0(cohort, "_", set_name, "_xCell_correlations_ALL.tsv")), sep = "\t")
    fwrite(cor_dt[FDR <= FDR_MAX],
           file.path(cohort_out, paste0(cohort, "_", set_name, "_xCell_correlations_FDR0.05.tsv")), sep = "\t")
    fwrite(cor_dt[FDR <= FDR_MAX & abs(rho) >= RHO_MIN],
           file.path(cohort_out, paste0(cohort, "_", set_name, "_xCell_correlations_FDR0.05_absRho0.3.tsv")), sep = "\t")
    
    # summary per score (percent of prioritized genes)
    sum_dt <- cor_dt[, .(
      n_genes_tested = uniqueN(gene),
      n_pos_strong   = sum(rho >=  RHO_MIN & FDR <= FDR_MAX),
      n_neg_strong   = sum(rho <= -RHO_MIN & FDR <= FDR_MAX)
    ), by = .(score_type)]
    
    sum_dt[, pct_pos := 100 * n_pos_strong / n_genes_tested]
    sum_dt[, pct_neg := 100 * n_neg_strong / n_genes_tested]
    
    fwrite(sum_dt, file.path(cohort_out, paste0(cohort, "_", set_name, "_xCell_directional_summary.tsv")), sep = "\t")
    
    # bar plot
    bar_dt <- rbind(
      sum_dt[, .(score_type, direction = "positive", pct = pct_pos)],
      sum_dt[, .(score_type, direction = "negative", pct = pct_neg)]
    )
    
    p_bar <- ggplot(bar_dt, aes(x = score_type, y = pct, fill = direction)) +
      geom_col() +
      labs(
        title = paste0(set_name, " prioritized HER2 features linked to admixture (", cohort, ")"),
        subtitle = paste0("Spearman |rho| ≥ ", RHO_MIN, ", BH-FDR ≤ ", FDR_MAX),
        x = "xCell score",
        y = "% of prioritized genes",
        fill = "Direction"
      ) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    
    ggsave(file.path(cohort_out, paste0(cohort, "_", set_name, "_xCell_barplot_directional.pdf")),
           p_bar, width = 7.6, height = 4.8)
    
    cor_dt
  }
  
  cor_p76 <- correlate_set(p76_genes, "P76")
  cor_p25 <- correlate_set(p25_genes, "P25")
  
  all_cor[[cohort]] <- rbindlist(list(cor_p76, cor_p25), fill = TRUE)
}

# pooled outputs
all_dt <- rbindlist(all_cor, fill = TRUE)
fwrite(all_dt, file.path(OUT_DIR, "ALLCOHORTS_P76P25_xCell_correlations_ALL.tsv"), sep = "\t")

# ===================== CROSS-COHORT SUMMARY PLOT =====================

# Defensive checks
all_dt <- all_dt[!is.na(cohort) & !is.na(set) & !is.na(gene) & !is.na(score_type)]
if (!nrow(all_dt)) stop("No pooled correlation results found; cannot build summary plot.")

# Summarize: within cohort x set x score_type
sum_all <- all_dt[, .(
  n_genes_tested = uniqueN(gene),
  n_pos_strong   = uniqueN(gene[rho >=  RHO_MIN & FDR <= FDR_MAX]),
  n_neg_strong   = uniqueN(gene[rho <= -RHO_MIN & FDR <= FDR_MAX])
), by = .(cohort, set, score_type)]

sum_all[, pct_pos := 100 * n_pos_strong / pmax(n_genes_tested, 1)]
sum_all[, pct_neg := 100 * n_neg_strong / pmax(n_genes_tested, 1)]

# Write summary table
fwrite(sum_all, file.path(OUT_DIR, "ALLCOHORTS_P76P25_xCell_directional_summary.tsv"), sep = "\t")

# Build plot data (stacked bars = positive + negative %)
bar_all <- rbind(
  sum_all[, .(cohort, set, score_type, direction = "positive", pct = pct_pos)],
  sum_all[, .(cohort, set, score_type, direction = "negative", pct = pct_neg)]
)

# Make ordering consistent (for slide readability)
bar_all[, cohort := factor(cohort, levels = c("Wolf","Robinson","Brueffer"))]
bar_all[, set := factor(set, levels = c("P76","P25"))]

p_summary <- ggplot(bar_all, aes(x = score_type, y = pct, fill = direction)) +
  geom_col() +
  facet_grid(set ~ cohort) +
  labs(
    title = "Prioritized HER2 feature sets show admixture-linked signal across cohorts",
    subtitle = paste0("Strong associations defined as |rho| ≥ ", RHO_MIN, " and BH-FDR ≤ ", FDR_MAX),
    x = "xCell score",
    y = "% of prioritized genes",
    fill = "Direction"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(OUT_DIR, "ALLCOHORTS_P76P25_xCell_summary_barplot.pdf"),
  plot = p_summary,
  width = 12.5, height = 7.2
)

cat("\nDONE. Outputs in: ", OUT_DIR, "\n", sep = "")
# =====================================================================
