# ================= phase1_Panels.R =================
# Makes a 3-panel volcano (A/B/C) from the cohort ALL_DE tables you already output.
# Output:
#   - Figure1_Volcano_3panel.pdf
#   - Figure1_Volcano_3panel.png

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  if (!requireNamespace("ggplot2", quietly = TRUE))    install.packages("ggplot2")
  if (!requireNamespace("patchwork", quietly = TRUE))  install.packages("patchwork")
  if (!requireNamespace("ggrepel", quietly = TRUE))    install.packages("ggrepel")
  if (!requireNamespace("ragg", quietly = TRUE))       install.packages("ragg")
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ragg)
})

# ---------------- USER CONFIG ----------------
wolf_all_fp <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/GSE194040_ISPY2_mRNA_BPsubtype/ALL_HER2pos_vs_HER2neg_DE.tsv"
rob_all_fp  <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/robinson et al/aim2_outputs/GSE199633_Robinson2025/ALL_HER2pos_vs_HER2neg_DE.tsv"
bru_all_fp  <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/brueffer et al/aim2_outputs/GSE81538/ALL_HER2pos_vs_HER2neg_DE.tsv"

out_dir <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/figures_phase1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Label genes to annotate (keep this SHORT or it gets messy)
LABEL_GENES <- c("ERBB2","GRB7","STARD3","MIEN1")

# Volcano thresholds (match what you describe in Results)
FDR_CUTOFF <- 0.05
LFC_CUTOFF <- 1.0

# ---------------- helpers ----------------
read_all_de <- function(fp) {
  stopifnot(file.exists(fp))
  dt <- fread(fp)
  # common columns from limma topTable:
  # logFC, P.Value, adj.P.Val, feature_id
  # Some of your files might use feature_id for gene symbol.
  if (!("feature_id" %in% names(dt))) {
    # try gene_symbol fallback
    if ("gene_symbol" %in% names(dt)) setnames(dt, "gene_symbol", "feature_id")
  }
  stopifnot(all(c("feature_id","logFC","adj.P.Val","P.Value") %in% names(dt)))
  dt[, feature_id := toupper(as.character(feature_id))]
  dt[, logFC := as.numeric(logFC)]
  dt[, adj.P.Val := as.numeric(adj.P.Val)]
  dt[, P.Value := as.numeric(P.Value)]
  dt
}

make_volcano_gg <- function(dt, title) {
  dt <- copy(dt)
  dt[, neglog10FDR := -log10(pmax(adj.P.Val, 1e-300))]
  dt[, sig := !is.na(adj.P.Val) & adj.P.Val < FDR_CUTOFF]
  dt[, sig_fc := sig & !is.na(logFC) & abs(logFC) >= LFC_CUTOFF]
  
  # Only label a few genes if present
  lab <- dt[feature_id %in% LABEL_GENES]
  
  ggplot(dt, aes(x = logFC, y = neglog10FDR)) +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = 2, linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = 2, linewidth = 0.3) +
    geom_point(aes(alpha = sig), size = 1.0) +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.95), guide = "none") +
    geom_point(data = dt[sig_fc == TRUE], size = 1.2) +
    ggrepel::geom_text_repel(
      data = lab,
      aes(label = feature_id),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.2,
      min.segment.length = 0
    ) +
    labs(
      title = title,
      x = "log2FC (HER2+ vs HER2-)",
      y = expression(-log[10]("FDR"))
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

# ---------------- build panels ----------------
wolf_dt <- read_all_de(wolf_all_fp)
rob_dt  <- read_all_de(rob_all_fp)
bru_dt  <- read_all_de(bru_all_fp)

p_wolf <- make_volcano_gg(wolf_dt, "Wolf (GSE194040)")
p_rob  <- make_volcano_gg(rob_dt,  "Robinson (GSE199633)")
p_bru  <- make_volcano_gg(bru_dt,  "Brueffer (GSE81538)")

panel <- (p_wolf | p_rob | p_bru) +
  plot_annotation(
    title = "HER2-associated transcriptional changes across cohorts",
    subtitle = sprintf("Volcano plots (FDR < %.2f; |log2FC| ≥ %.1f)", FDR_CUTOFF, LFC_CUTOFF),
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.99)
  )

# ---------------- save ----------------
png_out <- file.path(out_dir, "Figure1_Volcano_3panel.png")
pdf_out <- file.path(out_dir, "Figure1_Volcano_3panel.pdf")

ragg::agg_png(png_out, width = 3600, height = 1400, units = "px", res = 200, background = "white")
print(panel)
dev.off()

pdf(pdf_out, width = 18, height = 7, onefile = FALSE)
print(panel)
dev.off()

cat("Wrote:\n- ", png_out, "\n- ", pdf_out, "\n", sep = "")
# ==============================================================
# ================= Figure3_FDRonly_UpSetPlusTable.R =================
# BASE_DIR = your FDRonly folder (the one you just pasted)
# Panel A: UpSet (UP genes; FDR-only)
# Panel B: P76 / P25 summary table
#
# Writes into BASE_DIR:
#   Figure3_FDRonly_UpSetPlusTable.png
#   Figure3_FDRonly_UpSetPlusTable.pdf

suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
  library(ragg)
  library(png)
  library(readr)
})

# ---------------- CONFIG ----------------
BASE_DIR <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_upset_outputs/transcript_upset_dual_3cohorts/FDRonly"

UPSET_PNG <- file.path(BASE_DIR, "UpSet_UP_HER2pos_transcripts_FDR0.05.png")

# 🔴 Put your REAL Mac paths here (NOT /mnt/data)
P76_FILE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv"
P25_FILE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized/HER2_UP_overlap_25_FDR0.05_LFC1_ALL3_with_stats_Wolf_Robinson_Brueffer.csv"

OUT_PNG <- file.path(BASE_DIR, "Figure3_FDRonly_UpSetPlusTable.png")
OUT_PDF <- file.path(BASE_DIR, "Figure3_FDRonly_UpSetPlusTable.pdf")

# ---------------- helpers ----------------
stopifnot(dir.exists(BASE_DIR))
stopifnot(file.exists(UPSET_PNG))
stopifnot(file.exists(P76_FILE))
stopifnot(file.exists(P25_FILE))

add_panel_label <- function(g, label = "A") {
  grobTree(
    g,
    textGrob(label,
             x = unit(0.015, "npc"), y = unit(0.985, "npc"),
             just = c("left", "top"),
             gp = gpar(fontface = "bold", fontsize = 18))
  )
}

read_set_n <- function(fp) {
  x <- suppressWarnings(read_csv(fp, show_col_types = FALSE))
  if (!nrow(x)) return(0)
  if ("gene_symbol" %in% names(x)) return(length(unique(na.omit(x$gene_symbol))))
  length(unique(na.omit(x[[1]])))
}

# ---------------- Panel A ----------------
img <- png::readPNG(UPSET_PNG)
panelA <- rasterGrob(img, interpolate = TRUE)
panelA <- add_panel_label(panelA, "A")

# ---------------- Panel B ----------------
n_p76 <- read_set_n(P76_FILE)
n_p25 <- read_set_n(P25_FILE)

tab <- data.frame(
  Set = c("P76", "P25"),
  Definition = c(
    "Up in HER2+ vs HER2−, FDR < 0.05 in Wolf/Robinson/Brueffer (directionally consistent; 3-way overlap)",
    "P76 filtered to high-effect genes: FDR < 0.05 AND |log2FC| ≥ 1 in all 3 cohorts"
  ),
  `N genes` = c(n_p76, n_p25),
  check.names = FALSE
)

tg <- tableGrob(
  tab, rows = NULL,
  theme = ttheme_minimal(
    base_size = 14,
    core = list(fg_params = list(hjust = 0, x = 0.02)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.02))
  )
)

panelB <- add_panel_label(tg, "B")

# ---------------- combine ----------------
# Slightly more even balance than before: give Panel A more room, but not *all* the room.
fig <- arrangeGrob(
  panelA,
  panelB,
  ncol = 1,
  heights = c(3.0, 1.6)
)

# ---------------- write ----------------
ragg::agg_png(OUT_PNG, width = 2600, height = 2400, res = 200, background = "white")
grid.newpage(); grid.draw(fig)
dev.off()

pdf(OUT_PDF, width = 10.5, height = 10.0, onefile = FALSE)
grid.newpage(); grid.draw(fig)
dev.off()

message("Wrote: ", OUT_PNG)
message("Wrote: ", OUT_PDF)
# ===================== make_Figure6_proteomics_panel_1x3.R =====================

suppressPackageStartupMessages({
  library(pdftools)
  library(magick)
  library(grid)
})

# ---- INPUT PDF FILES ----
pdf_A <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/raj kumar et al/aim2_outputs/RajKumar_PROT_GLOBAL/volcano.pdf"
pdf_B <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/mertins et al/aim2_outputs/Mertins_CPTAC_BRCA_PROT_GLOBAL/volcano.pdf"
pdf_C <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/krug et al/aim2_outputs/Krug_CPTAC_BRCA_PROT_A_ClinicalOnly/volcano.pdf"

# ---- OUTPUT FILES ----
out_dir <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA"
out_png <- file.path(out_dir, "Figure6_ProteomicsVolcano_ABC.png")
out_pdf <- file.path(out_dir, "Figure6_ProteomicsVolcano_ABC.pdf")

# ---- Convert first page of PDF to PNG ----
pdf1_to_png <- function(pdf_path, out_prefix, dpi = 300) {
  if (!file.exists(pdf_path)) stop("PDF not found:\n  ", pdf_path)
  
  out_files <- pdftools::pdf_convert(
    pdf = pdf_path,
    format = "png",
    pages = 1,
    dpi = dpi,
    filenames = paste0(out_prefix, "_page")
  )
  
  if (length(out_files) < 1 || !file.exists(out_files[1])) {
    stop("PDF conversion failed for:\n  ", pdf_path)
  }
  out_files[1]
}

# ---- Make label as transparent PNG (avoids ImageMagick font bugs) ----
make_label_png <- function(label, px = 140, fontsize = 90) {
  tf <- tempfile(fileext = ".png")
  png(tf, width = px, height = px, bg = "transparent")
  grid.newpage()
  grid.text(
    label,
    x = unit(0.1, "npc"),
    y = unit(0.9, "npc"),
    just = c("left", "top"),
    gp = gpar(col = "black", fontsize = fontsize, fontface = "bold")
  )
  dev.off()
  tf
}

add_label <- function(im, label = "A", border_px = 40, x_off = 20, y_off = 20) {
  im2 <- image_border(im, color = "white", geometry = paste0(border_px, "x", border_px))
  lab_file <- make_label_png(label)
  lab_img  <- image_read(lab_file)
  
  im3 <- image_composite(
    im2, lab_img,
    operator = "over",
    gravity = "northwest",
    offset = paste0("+", x_off, "+", y_off)
  )
  
  unlink(lab_file)
  im3
}

# ---- Convert PDFs to PNGs ----
png_A <- pdf1_to_png(pdf_A, file.path(out_dir, "tmp_Fig6_A"))
png_B <- pdf1_to_png(pdf_B, file.path(out_dir, "tmp_Fig6_B"))
png_C <- pdf1_to_png(pdf_C, file.path(out_dir, "tmp_Fig6_C"))

# ---- Read PNGs ----
A <- image_read(png_A)
B <- image_read(png_B)
C <- image_read(png_C)

# ---- Normalize sizes (controls final figure size) ----
target_w <- 1200   # increase to 1400–1600 if you want bigger
A <- image_resize(A, paste0(target_w, "x"))
B <- image_resize(B, paste0(target_w, "x"))
C <- image_resize(C, paste0(target_w, "x"))

# ---- Add A/B/C labels ----
A <- add_label(A, "A")
B <- add_label(B, "B")
C <- add_label(C, "C")

# ---- Assemble single-row panel ----
row_ABC <- image_append(c(A, B, C), stack = FALSE)

final <- image_border(row_ABC, color = "white", geometry = "40x40")

# ---- Write outputs ----
image_write(final, path = out_png, format = "png", density = 300)
image_write(final, path = out_pdf, format = "pdf", density = 300)

cat("✅ Saved:\n", out_png, "\n", out_pdf, "\n", sep = "")
# =============================================================================
# ================= Prioritized (P76 / P25) RNA–Protein Concordance =================
# FINAL (EDITED): cleaner PATH config + auto 2-panel Figure (P76 + P25) per mode
# - Uses RNA consensus (median logFC across Wolf/Robinson/Brueffer)
# - Tests protein support (RajKumar global proteome)
# - Writes per-set tables + per-set scatter + a combined 2-panel figure (A,B)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(patchwork)   # <- for 2-panel figure
})

# ---------------- CONFIG ----------------
FDR_MAX <- 0.05
FC_MIN  <- 1

# ---- BASES (EDIT THESE ONCE, EVERYTHING ELSE DERIVES FROM THEM) ----
THESIS_BASE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA"

OUT_BASE <- file.path(THESIS_BASE, "aim2_concordance_outputs", "Prioritized_RNA_PROT_FINAL")
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

# ---------------- INPUTS ----------------
P76_FILE <- file.path(
  THESIS_BASE, "aim2_prioritized",
  "HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv"
)

P25_FILE <- file.path(
  THESIS_BASE, "aim2_prioritized",
  "HER2_UP_overlap_25_FDR0.05_LFC1_ALL3_with_stats_Wolf_Robinson_Brueffer.csv"
)

WOLF_RNA <- file.path(
  THESIS_BASE, "wolf et al", "aim2_outputs",
  "GSE194040_ISPY2_mRNA_BPsubtype", "ALL_HER2pos_vs_HER2neg_DE.tsv"
)

ROBINSON_RNA <- file.path(
  THESIS_BASE, "robinson et al", "aim2_outputs",
  "GSE199633_Robinson2025", "ALL_HER2pos_vs_HER2neg_DE.tsv"
)

BRUEFFER_RNA <- file.path(
  THESIS_BASE, "brueffer et al", "aim2_outputs",
  "GSE81538", "ALL_HER2pos_vs_HER2neg_DE.tsv"
)

PROT_FILE <- file.path(
  THESIS_BASE, "raj kumar et al", "aim2_outputs",
  "RajKumar_PROT_GLOBAL", "DE_full.tsv"
)

# ---------------- HELPERS ----------------
clean_gene <- function(x) toupper(str_trim(gsub("\\.\\d+$", "", as.character(x))))
sign_dir   <- function(x) ifelse(x > 0, "UP", ifelse(x < 0, "DOWN", NA))

canon <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

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
    cat("File:", path, "\nColumns seen:\n")
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
    cat("File:", path, "\nColumns seen:\n")
    print(names(dt))
    stop("Cannot detect columns in RajKumar PROT DE.")
  }
  
  out <- dt[, .(
    gene      = clean_gene(get(gene_col)),
    prot_logFC = as.numeric(get(fc_col)),
    prot_fdr   = as.numeric(get(fdr_col))
  )]
  out <- out[!is.na(gene) & gene != "" & !is.na(prot_logFC) & !is.na(prot_fdr)]
  out <- out[prot_fdr <= FDR_MAX]
  out[, prot_dir := sign_dir(prot_logFC)]
  out <- out[!is.na(prot_dir)]
  unique(out)
}

# ---------------- PLOTTING (returns a ggplot OR NULL) ----------------
make_scatter <- function(m, set_name, mode) {
  if (!nrow(m)) return(NULL)
  
  label_genes <- intersect(
    c("ERBB2","GRB7","STARD3","PGAP3","MIEN1","CDK12","CCNK","ERBB4"),
    m$gene
  )
  
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
  
  p
}

# ---------------- MAIN: analyze one set ----------------
analyze_set <- function(set_name, genes, rna_all, rna_cons, prot, outdir, mode) {
  
  rna_set  <- rna_cons[gene %in% genes]
  prot_set <- prot[gene %in% genes]
  
  m <- merge(rna_set, prot_set, by = "gene")
  m[, concordant := (rna_dir == prot_dir)]
  
  # tables
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
  
  # RNA support count among overlapped genes
  if (nrow(m) > 0) {
    rna_support <- rna_all[gene %in% m$gene, .(rna_support_n = uniqueN(cohort)), by = gene]
    m2 <- merge(m, rna_support, by = "gene", all.x = TRUE)
    fwrite(m2[order(-rna_support_n, -abs(prot_logFC))],
           file.path(outdir, paste0(set_name, "_overlap_with_RNA_support.tsv")),
           sep = "\t")
  }
  
  # plot (save single)
  p <- make_scatter(m, set_name, mode)
  if (!is.null(p)) {
    ggsave(file.path(outdir, paste0(set_name, "_scatter.png")),
           p, width = 10, height = 6, dpi = 300)
    ggsave(file.path(outdir, paste0(set_name, "_scatter.pdf")),
           p, width = 10, height = 6)
  }
  
  list(summary = summary, plot = p)
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
    rna_logFC   = median(logFC, na.rm = TRUE),
    rna_fdr_min = min(fdr, na.rm = TRUE),
    n_cohorts   = uniqueN(cohort)
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
  a1 <- analyze_set("P76", P76, rna_all, rna_cons, prot, outdir, mode)
  a2 <- analyze_set("P25", P25, rna_all, rna_cons, prot, outdir, mode)
  
  # combined summary table
  fwrite(rbindlist(list(a1$summary, a2$summary), fill = TRUE),
         file.path(outdir, "SUMMARY_P76_P25.tsv"),
         sep = "\t")
  
  # ---- NEW: 2-panel figure (A,B) ----
  # If either plot is NULL (no overlap), we still write something sensible.
  if (!is.null(a1$plot) && !is.null(a2$plot)) {
    fig2 <- (a1$plot + a2$plot) +
      plot_layout(ncol = 2) +
      plot_annotation(tag_levels = "A")  # adds A, B automatically
    
    ggsave(file.path(outdir, "P76_P25_scatter_2panel.png"),
           fig2, width = 18, height = 6, dpi = 300)
    ggsave(file.path(outdir, "P76_P25_scatter_2panel.pdf"),
           fig2, width = 18, height = 6)
  } else {
    message("NOTE: Could not make 2-panel figure for ", mode,
            " (one of P76/P25 had 0 overlap).")
  }
  
  cat("\nDONE:", mode, "->", outdir, "\n")
}

# ---------------- RUN BOTH MODES ----------------
run_mode("FDR_only")
run_mode("FDR_FC")

cat("\nALL DONE -> ", OUT_BASE, "\n", sep = "")
