# TranscriptUpSetPlot.R — FINAL (3 MAIN TRANSCRIPT COHORTS ONLY)
# DUAL OUTPUTS: FDR-only + FDR+FC
#  + UpSet/Venn
#  + Cross-cohort concordance (ALL thresholded genes) with highlighted all-3 shared genes
#  + Concordance 3-panel tagged A/B/C

suppressPackageStartupMessages({
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  if (!requireNamespace("data.table", quietly = TRUE))        install.packages("data.table")
  if (!requireNamespace("ragg", quietly = TRUE))             install.packages("ragg")
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))   install.packages("ComplexHeatmap")
  if (!requireNamespace("circlize", quietly = TRUE))         install.packages("circlize")
  if (!requireNamespace("VennDiagram", quietly = TRUE))      install.packages("VennDiagram")
  if (!requireNamespace("ggplot2", quietly = TRUE))          install.packages("ggplot2")
  if (!requireNamespace("patchwork", quietly = TRUE))        install.packages("patchwork")
  if (!requireNamespace("ggrepel", quietly = TRUE))          install.packages("ggrepel")
  
  library(data.table)
})

# ---------------- CONFIG ----------------
FDR_MAX <- 0.05

MODES <- list(
  FDRonly = list(tag = "FDR0.05",      lfc_min = 0.00, prefer_lfc_file = FALSE),
  FDR_FC  = list(tag = "FDR0.05_LFC1", lfc_min = 1.00, prefer_lfc_file = TRUE)
)

base_out <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_upset_outputs"
out_dir  <- file.path(base_out, "transcript_upset_dual_3cohorts")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 3 cohorts only
roots_main <- c(
  "GSE194040_ISPY2_mRNA_BPsubtype" = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/GSE194040_ISPY2_mRNA_BPsubtype",
  "GSE199633_Robinson2025"         = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/robinson et al/aim2_outputs/GSE199633_Robinson2025",
  "GSE81538"                       = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/brueffer et al/aim2_outputs/GSE81538"
)

BIO_ORDER <- c("GSE 194040", "GSE 199633", "GSE 81538")

# ---------------- helpers ----------------
pick_de_file <- function(pathlike, prefer_lfc_file = FALSE) {
  if (file.exists(pathlike) && !dir.exists(pathlike)) return(pathlike)
  if (!dir.exists(pathlike)) return(NA_character_)
  files <- list.files(pathlike, full.names = TRUE, recursive = TRUE)
  
  pref_fdr <- grep("FILT_HER2pos_vs_HER2neg_DE_FDR0\\.05\\.(tsv|csv)$",
                   files, ignore.case = TRUE, value = TRUE)
  pref_lfc <- grep("FILT_HER2pos_vs_HER2neg_DE_FDR0\\.05_LFC1\\.(tsv|csv)$",
                   files, ignore.case = TRUE, value = TRUE)
  
  if (prefer_lfc_file) {
    if (length(pref_lfc)) return(pref_lfc[1])
    if (length(pref_fdr)) return(pref_fdr[1])
  } else {
    if (length(pref_fdr)) return(pref_fdr[1])
    if (length(pref_lfc)) return(pref_lfc[1])
  }
  
  cand <- grep("HER2pos_vs_HER2neg.*(DE|DEP).*\\.(tsv|csv)$",
               files, ignore.case = TRUE, value = TRUE)
  if (length(cand)) return(cand[1])
  
  NA_character_
}

read_tab <- function(p) {
  if (grepl("\\.csv$", p, ignore.case = TRUE)) {
    read.csv(p, check.names = FALSE)
  } else {
    read.table(p, header = TRUE, sep = "\t", quote = "",
               comment.char = "", check.names = FALSE)
  }
}

pick_col <- function(cols, cands, allow_prefix = TRUE) {
  hit <- intersect(cols, cands); if (length(hit)) return(hit[1])
  low <- tolower(cols)
  for (c in cands) {
    i <- which(low == tolower(c)); if (length(i)) return(cols[i[1]])
  }
  if (allow_prefix) for (c in cands) {
    i <- which(startsWith(low, tolower(c))); if (length(i)) return(cols[i[1]])
  }
  NA_character_
}

clean_syms <- function(v) {
  v <- as.character(v)
  v <- gsub("^\\s+|\\s+$", "", v)
  v <- gsub("\\.\\d+$", "", v)
  toupper(v)
}

pretty_names <- function(x) {
  x <- gsub("_.*$", "", x)
  x <- gsub("^GSE", "GSE ", x)
  x
}

# ---- strict overlap sets (directionally consistent by definition) ----
get_overlap_up_down <- function(sets_up, sets_down) {
  sets_up   <- lapply(sets_up,   function(x) unique(toupper(as.character(x))))
  sets_down <- lapply(sets_down, function(x) unique(toupper(as.character(x))))
  
  up_all3   <- Reduce(intersect, sets_up)
  down_all3 <- Reduce(intersect, sets_down)
  
  conflict <- intersect(up_all3, down_all3)
  if (length(conflict)) {
    up_all3   <- setdiff(up_all3, conflict)
    down_all3 <- setdiff(down_all3, conflict)
  }
  list(
    up = sort(unique(up_all3)),
    down = sort(unique(down_all3)),
    conflict = sort(unique(conflict))
  )
}

# ---------------- UpSet + Venn ----------------
make_upset <- function(sets_list, title_base, file_stub, out_dir_mode) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap); library(circlize); library(grid); library(ragg)
  })
  cm <- ComplexHeatmap::make_comb_mat(sets_list)
  if (length(ComplexHeatmap::comb_size(cm)) == 0) return(invisible(NULL))
  
  top_idx <- order(ComplexHeatmap::comb_size(cm), decreasing = TRUE)
  top_idx <- top_idx[seq_len(min(30, length(top_idx)))]
  
  ht <- ComplexHeatmap::UpSet(
    cm[top_idx],
    pt_size = unit(3, "mm"), lwd = 1.2,
    top_annotation = ComplexHeatmap::upset_top_annotation(
      cm[top_idx], add_numbers = TRUE, gp = gpar(fill = "gray25")
    ),
    set_order = seq_along(ComplexHeatmap::set_name(cm))
  )
  
  ragg::agg_png(file.path(out_dir_mode, paste0("UpSet_", file_stub, ".png")),
                width = 2600, height = 1400, units = "px", res = 200, background = "white")
  grid::grid.newpage()
  grid::grid.text(title_base, y = unit(1, "npc") - unit(8, "pt"))
  ComplexHeatmap::draw(ht, padding = unit(c(10,10,10,10), "pt"))
  dev.off()
  
  pdf(file.path(out_dir_mode, paste0("UpSet_", file_stub, ".pdf")),
      width = 17, height = 10, onefile = FALSE)
  grid::grid.newpage()
  grid::grid.text(title_base, y = unit(1, "npc") - unit(8, "pt"))
  ComplexHeatmap::draw(ht, padding = unit(c(10,10,10,10), "pt"))
  dev.off()
}

make_venn3 <- function(sets_list, title, stub, out_dir_mode) {
  if (length(sets_list) != 3) return(invisible(NULL))
  suppressPackageStartupMessages({ library(VennDiagram); library(grid); library(ragg) })
  
  nms <- names(sets_list)
  n12  <- length(intersect(sets_list[[1]], sets_list[[2]]))
  n13  <- length(intersect(sets_list[[1]], sets_list[[3]]))
  n23  <- length(intersect(sets_list[[2]], sets_list[[3]]))
  n123 <- length(Reduce(intersect, sets_list))
  
  ragg::agg_png(file.path(out_dir_mode, paste0("Venn_", stub, ".png")),
                width = 1800, height = 1400, units = "px", res = 200, background = "white")
  grid::grid.newpage()
  grid::grid.text(title, y = unit(1, "npc") - unit(10, "pt"))
  vp <- viewport(y = 0.48, height = 0.9); pushViewport(vp)
  VennDiagram::draw.triple.venn(
    area1 = length(sets_list[[1]]),
    area2 = length(sets_list[[2]]),
    area3 = length(sets_list[[3]]),
    n12 = n12, n13 = n13, n23 = n23, n123 = n123,
    category = nms,
    fill = c("grey85", "grey70", "grey55"),
    alpha = rep(0.7, 3),
    cex = 1.2, cat.cex = 1.2, cat.pos = c(-20, 20, 0)
  )
  popViewport()
  dev.off()
  
  pdf(file.path(out_dir_mode, paste0("Venn_", stub, ".pdf")), width = 9, height = 7)
  grid::grid.newpage()
  grid::grid.text(title, y = unit(1, "npc") - unit(10, "pt"))
  vp <- viewport(y = 0.48, height = 0.9); pushViewport(vp)
  VennDiagram::draw.triple.venn(
    area1 = length(sets_list[[1]]),
    area2 = length(sets_list[[2]]),
    area3 = length(sets_list[[3]]),
    n12 = n12, n13 = n13, n23 = n23, n123 = n123,
    category = nms,
    fill = c("grey85", "grey70", "grey55"),
    alpha = rep(0.7, 3),
    cex = 1.2, cat.cex = 1.2, cat.pos = c(-20, 20, 0)
  )
  popViewport()
  dev.off()
}

# ---------------- cross-cohort concordance ----------------
make_pair_concordance_allgenes <- function(de_tables, thr_sets_union, highlight_genes,
                                           tag, x_study, y_study,
                                           label_genes = character(0)) {
  suppressPackageStartupMessages({ library(ggplot2) })
  
  if (!(x_study %in% names(de_tables)) || !(y_study %in% names(de_tables))) return(NULL)
  
  x <- as.data.table(de_tables[[x_study]])
  y <- as.data.table(de_tables[[y_study]])
  
  # only genes that passed thresholds (UP ∪ DOWN) in EACH cohort for this mode
  x_keep <- intersect(x$gene_symbol, thr_sets_union[[x_study]])
  y_keep <- intersect(y$gene_symbol, thr_sets_union[[y_study]])
  keep   <- intersect(x_keep, y_keep)
  
  if (length(keep) < 10) return(NULL)
  
  x <- x[gene_symbol %in% keep]
  y <- y[gene_symbol %in% keep]
  
  m <- merge(x[, .(gene_symbol, x_log2FC = log2FC)],
             y[, .(gene_symbol, y_log2FC = log2FC)],
             by = "gene_symbol", all = FALSE)
  
  if (nrow(m) < 10) return(NULL)
  
  m[, is_shared_all3 := gene_symbol %in% highlight_genes]
  
  sp <- suppressWarnings(cor(m$x_log2FC, m$y_log2FC, method = "spearman", use = "complete.obs"))
  pr <- suppressWarnings(cor(m$x_log2FC, m$y_log2FC, method = "pearson",  use = "complete.obs"))
  
  ttl <- sprintf("%s vs %s", x_study, y_study)
  sub <- sprintf("Spearman ρ=%.2f | Pearson r=%.2f | n=%d genes", sp, pr, nrow(m))
  
  p <- ggplot(m, aes(x = x_log2FC, y = y_log2FC)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(aes(alpha = is_shared_all3), size = 1.4) +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.95), guide = "none") +
    geom_point(data = m[is_shared_all3 == TRUE], size = 2.1) +
    labs(
      title = ttl,
      subtitle = sub,
      x = sprintf("log2FC (HER2+ vs HER2-) in %s", x_study),
      y = sprintf("log2FC (HER2+ vs HER2-) in %s", y_study)
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  
  if (length(label_genes) > 0) {
    suppressPackageStartupMessages(library(ggrepel))
    lab <- m[gene_symbol %in% toupper(label_genes)]
    if (nrow(lab) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = lab,
        aes(label = gene_symbol),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.3,
        point.padding = 0.2,
        min.segment.length = 0
      )
    }
  }
  
  p
}

make_concordance_3panel <- function(de_tables, thr_sets_union, highlight_genes,
                                    out_dir_mode, tag,
                                    label_genes = c("ERBB2","GRB7","STARD3","MIEN1")) {
  suppressPackageStartupMessages({ library(patchwork); library(ragg) })
  
  pairs <- list(
    c("GSE 194040", "GSE 199633"),
    c("GSE 194040", "GSE 81538"),
    c("GSE 199633", "GSE 81538")
  )
  
  plots <- list()
  for (pr in pairs) {
    plots[[paste(pr, collapse = "_vs_")]] <-
      make_pair_concordance_allgenes(
        de_tables = de_tables,
        thr_sets_union = thr_sets_union,
        highlight_genes = highlight_genes,
        tag = tag,
        x_study = pr[1],
        y_study = pr[2],
        label_genes = label_genes
      )
  }
  
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) return(invisible(NULL))
  
  fig <- wrap_plots(plots, nrow = 1) +
    plot_annotation(
      title = sprintf("Cross-cohort concordance of HER2-associated effect sizes (%s)", tag),
      subtitle = "All thresholded genes shown; genes shared across all 3 cohorts are highlighted",
      tag_levels = "A"
    ) &
    theme(
      plot.tag = element_text(face = "bold", size = 16),
      plot.tag.position = c(0.01, 0.99)
    )
  
  stub <- paste0("CrossCohortConcordance_3panel_", tag)
  
  png_out <- file.path(out_dir_mode, paste0(stub, ".png"))
  ragg::agg_png(png_out, width = 3600, height = 1400, units = "px", res = 200, background = "white")
  print(fig); dev.off()
  
  pdf_out <- file.path(out_dir_mode, paste0(stub, ".pdf"))
  pdf(pdf_out, width = 18, height = 7, onefile = FALSE)
  print(fig); dev.off()
  
  cat("Wrote concordance figure: ", png_out, "\n", sep = "")
}

# ---------------- RUN BOTH MODES ----------------
for (mode_name in names(MODES)) {
  mode <- MODES[[mode_name]]
  LFC_MIN <- mode$lfc_min
  
  out_dir_mode <- file.path(out_dir, mode_name)
  dir.create(out_dir_mode, showWarnings = FALSE, recursive = TRUE)
  
  de_files <- vapply(
    roots_main,
    pick_de_file,
    FUN.VALUE = character(1),
    prefer_lfc_file = mode$prefer_lfc_file
  )
  if (anyNA(de_files)) {
    bad <- names(de_files)[is.na(de_files)]
    stop("No DE file found in:\n- ", paste(roots_main[bad], collapse = "\n- "))
  }
  
  audit_rows <- list()
  sets_up <- list()
  sets_down <- list()
  de_tables <- list()
  
  for (lbl in names(de_files)) {
    p  <- de_files[[lbl]]
    dt <- read_tab(p)
    
    g <- pick_col(names(dt), c("gene_symbol","SYMBOL","symbol","Gene","GeneSymbol","gene","feature_id"))
    q <- pick_col(names(dt), c("adj.P.Val","padj","FDR","qval","adj_pval","FDR.BH"))
    l <- pick_col(names(dt), c("logFC","log2FoldChange","log2FC","log2FC_HER2pos_vs_HER2neg"))
    
    if (any(is.na(c(g,q,l)))) {
      audit_rows[[lbl]] <- data.frame(study=lbl, file=p, n_total=nrow(dt),
                                      n_sig=NA, n_up=NA, n_down=NA, note="missing_cols")
      next
    }
    
    fdr <- suppressWarnings(as.numeric(dt[[q]]))
    lfc <- suppressWarnings(as.numeric(dt[[l]]))
    sym <- clean_syms(dt[[g]])
    
    ok <- which(!is.na(fdr) & !is.na(lfc) & !is.na(sym) &
                  fdr <= FDR_MAX & abs(lfc) >= LFC_MIN)
    
    up_syms   <- unique(sym[ok][lfc[ok] > 0])
    down_syms <- unique(sym[ok][lfc[ok] < 0])
    
    sets_up[[lbl]]   <- up_syms
    sets_down[[lbl]] <- down_syms
    
    audit_rows[[lbl]] <- data.frame(
      study=lbl, file=p, n_total=nrow(dt),
      n_sig=length(ok), n_up=length(up_syms), n_down=length(down_syms),
      note=""
    )
    
    study_nm <- pretty_names(lbl)
    de_tables[[study_nm]] <- data.table(gene_symbol = sym, log2FC = lfc, FDR = fdr)
  }
  
  audit_df <- data.table::rbindlist(audit_rows, use.names = TRUE, fill = TRUE)
  write.csv(audit_df,
            file.path(out_dir_mode, paste0("DE_audit_counts_transcript_", mode$tag, ".csv")),
            row.names = FALSE)
  
  names(sets_up)   <- pretty_names(names(sets_up))
  names(sets_down) <- pretty_names(names(sets_down))
  
  sets_up   <- sets_up[BIO_ORDER]
  sets_down <- sets_down[BIO_ORDER]
  
  # per-cohort lists
  for (nm in names(sets_up)) {
    writeLines(sort(sets_up[[nm]]),   file.path(out_dir_mode, paste0(nm, "_UP_", mode$tag, ".txt")))
    writeLines(sort(sets_down[[nm]]), file.path(out_dir_mode, paste0(nm, "_DOWN_", mode$tag, ".txt")))
  }
  
  # overlap UP/DOWN (directionally consistent across all 3)
  ov <- get_overlap_up_down(sets_up, sets_down)
  up_all3   <- ov$up
  down_all3 <- ov$down
  
  fwrite(data.table(gene_symbol = up_all3, direction = "UP"),
         file.path(out_dir_mode, paste0("HER2_Prioritized_UP_Overlap_Phase1_", mode$tag, ".csv")))
  fwrite(data.table(gene_symbol = down_all3, direction = "DOWN"),
         file.path(out_dir_mode, paste0("HER2_Prioritized_DOWN_Overlap_Phase1_", mode$tag, ".csv")))
  
  summary_dt <- data.table(
    mode = mode_name,
    tag = mode$tag,
    n_up_all3 = length(up_all3),
    n_down_all3 = length(down_all3),
    n_total = length(up_all3) + length(down_all3),
    n_conflicts_removed = length(ov$conflict)
  )
  fwrite(summary_dt, file.path(out_dir_mode, paste0("HER2_Prioritized_Overlap_", mode$tag, "_SUMMARY.csv")))
  
  cat("UP overlap (all 3): ", length(up_all3), "\n", sep = "")
  cat("DOWN overlap (all 3): ", length(down_all3), "\n", sep = "")
  
  # thresholded gene sets per cohort for this mode (UP ∪ DOWN)
  thr_sets_union <- lapply(names(sets_up), function(nm) unique(c(sets_up[[nm]], sets_down[[nm]])))
  names(thr_sets_union) <- names(sets_up)
  
  # highlight genes = (all-3 UP overlap) ∪ (all-3 DOWN overlap)
  highlight_genes <- unique(c(up_all3, down_all3))
  
  # concordance figure (3 panels) — NOW TAGGED A/B/C
  make_concordance_3panel(
    de_tables = de_tables,
    thr_sets_union = thr_sets_union,
    highlight_genes = highlight_genes,
    out_dir_mode = out_dir_mode,
    tag = mode$tag,
    label_genes = c("ERBB2","GRB7","STARD3","MIEN1")
  )
  
  # UpSet + Venn
  suffix <- mode$tag
  
  make_upset(
    sets_up,
    if (LFC_MIN > 0)
      sprintf("Up in HER2+ (transcripts; FDR≤%.2f & |log2FC|≥%.1f)", FDR_MAX, LFC_MIN)
    else
      sprintf("Up in HER2+ (transcripts; FDR≤%.2f)", FDR_MAX),
    paste0("UP_HER2pos_transcripts_", suffix),
    out_dir_mode
  )
  
  make_upset(
    sets_down,
    if (LFC_MIN > 0)
      sprintf("Down in HER2+ (transcripts; FDR≤%.2f & |log2FC|≥%.1f)", FDR_MAX, LFC_MIN)
    else
      sprintf("Down in HER2+ (transcripts; FDR≤%.2f)", FDR_MAX),
    paste0("DOWN_HER2pos_transcripts_", suffix),
    out_dir_mode
  )
  
  make_venn3(
    sets_up,
    if (LFC_MIN > 0)
      sprintf("3-way Venn — Up in HER2+ (FDR≤%.2f & |log2FC|≥%.1f)", FDR_MAX, LFC_MIN)
    else
      sprintf("3-way Venn — Up in HER2+ (FDR≤%.2f)", FDR_MAX),
    paste0("UP_HER2pos_", suffix),
    out_dir_mode
  )
  
  make_venn3(
    sets_down,
    if (LFC_MIN > 0)
      sprintf("3-way Venn — Down in HER2+ (FDR≤%.2f & |log2FC|≥%.1f)", FDR_MAX, LFC_MIN)
    else
      sprintf("3-way Venn — Down in HER2+ (FDR≤%.2f)", FDR_MAX),
    paste0("DOWN_HER2pos_", suffix),
    out_dir_mode
  )
  
  cat("Saved outputs → ", out_dir_mode, "\n", sep = "")
}

cat("DONE. Dual outputs saved under: ", out_dir, "\n", sep = "")
