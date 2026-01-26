# Aim2UpSetPlot_PROTEOME_overlap_FINAL.R — FINAL
# (3 PROTEOMICS COHORTS ONLY: Rajkumar / Mertins / Krug)
# Outputs:
#   - UpSet (UP + DOWN) for 2 regimes: FDR-only and FDR+FC
#   - 3-way Venn (UP + DOWN) for both regimes
#   - Audit CSV of parsed files and set sizes

suppressPackageStartupMessages({
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  if (!requireNamespace("ragg", quietly = TRUE))        install.packages("ragg")
  library(data.table)
})

# ---------------- CONFIG ----------------
FDR_MAX <- 0.05
LFC_MIN <- 1.00   # |log2FC| >= 1 matches your transcript strict regime

OUT_BASE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_upset_outputs"
OUT_DIR  <- file.path(OUT_BASE, "proteome_overlap_upset")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Restrict to the 3 proteomic cohorts you actually describe in Methods
roots <- c(
  RajKumar = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/raj kumar et al/aim2_outputs",
  Mertins  = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/mertins et al/aim2_outputs",
  Krug     = "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/krug et al/aim2_outputs"
)

# ---------------- HELPERS ----------------
is_proteomeish_path <- function(p) {
  grepl("PROT|PROTEOM|protein|global_prot|proteome", p, ignore.case = TRUE)
}

read_tab <- function(p) {
  if (grepl("\\.csv$", p, ignore.case = TRUE)) {
    read.csv(p, check.names = FALSE)
  } else {
    read.table(p, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
  }
}

pick_col <- function(cols, cands, allow_prefix = TRUE) {
  hit <- intersect(cols, cands); if (length(hit)) return(hit[1])
  low <- tolower(cols)
  for (c in cands) {
    i <- which(low == tolower(c)); if (length(i)) return(cols[i[1]])
  }
  if (allow_prefix) {
    for (c in cands) {
      i <- which(startsWith(low, tolower(c)))
      if (length(i)) return(cols[i[1]])
    }
  }
  NA_character_
}

clean_syms_prot <- function(v) {
  v <- as.character(v)
  v <- gsub("^\\s+|\\s+$", "", v)
  v <- gsub("\\.\\d+$", "", v)                          # Ensembl-like version tails
  v <- gsub("\\.total$", "", v, ignore.case = TRUE)     # RPPA-ish suffix
  v <- gsub("_TOTAL$", "", v, ignore.case = TRUE)
  v <- gsub("\\.pY\\d+$", "", v, ignore.case = TRUE)    # phospho-ish tags
  v <- gsub("[- ]", "_", v)
  toupper(v)
}

looks_like_dep_table <- function(df) {
  cn <- tolower(names(df))
  has_fc  <- any(c("logfc","log2foldchange","log2fc","log2fc_her2pos_vs_her2neg",
                   "log2fc.her2pos_vs_her2neg","lfc") %in% cn)
  has_fdr <- any(c("adj.p.val","padj","fdr","qval","adj_pval","fdr.bh","adj.pval") %in% cn)
  has_id  <- any(c("gene_symbol","symbol","gene","genesymbol","feature_id",
                   "protein","protein_id","hgnc_symbol") %in% cn)
  has_fc && has_fdr && has_id
}

make_sets_from_de <- function(df, fdr_max, lfc_min = NULL) {
  # detect columns
  g <- pick_col(names(df), c("gene_symbol","SYMBOL","symbol","Gene","GeneSymbol","gene",
                             "feature_id","Protein","protein","protein_id","HGNC_symbol","hgnc_symbol"))
  q <- pick_col(names(df), c("adj.P.Val","adj.Pval","adj.p.val","padj","FDR","qval","adj_pval","FDR.BH"))
  l <- pick_col(names(df), c("logFC","log2FoldChange","log2FC","log2FC_HER2pos_vs_HER2neg",
                             "log2FC.her2pos_vs_her2neg","lfc"))
  if (any(is.na(c(g, q, l)))) return(list(up = character(0), down = character(0), n_keep = 0))
  
  fdr  <- suppressWarnings(as.numeric(df[[q]]))
  lfc  <- suppressWarnings(as.numeric(df[[l]]))
  sym  <- clean_syms_prot(df[[g]])
  
  keep <- which(!is.na(fdr) & !is.na(lfc) & !is.na(sym) & fdr <= fdr_max)
  if (!is.null(lfc_min)) {
    keep <- keep[abs(lfc[keep]) >= lfc_min]
  }
  
  up   <- unique(sym[keep][lfc[keep] > 0])
  down <- unique(sym[keep][lfc[keep] < 0])
  
  list(up = up, down = down, n_keep = length(keep))
}

# ---------------- FIND CANDIDATE FILES ----------------
cand_files <- character(0)
for (r in unname(roots)) {
  if (!dir.exists(r)) next
  f <- list.files(r, pattern = "\\.(tsv|csv)$", full.names = TRUE, recursive = TRUE)
  f <- f[is_proteomeish_path(f) | is_proteomeish_path(dirname(f))]
  cand_files <- c(cand_files, f)
}
cand_files <- unique(cand_files)

# prefer explicit HER2pos_vs_HER2neg protein tables first
priority <- c(
  grep("DEP.*HER2pos_vs_HER2neg\\.(tsv|csv)$", cand_files, value = TRUE, ignore.case = TRUE),
  grep("DE.*HER2pos_vs_HER2neg\\.(tsv|csv)$",  cand_files, value = TRUE, ignore.case = TRUE),
  setdiff(cand_files, c(
    grep("DEP.*HER2pos_vs_HER2neg\\.(tsv|csv)$", cand_files, value = TRUE, ignore.case = TRUE),
    grep("DE.*HER2pos_vs_HER2neg\\.(tsv|csv)$",  cand_files, value = TRUE, ignore.case = TRUE)
  ))
)
cand_files <- unique(priority)

# label cohort based on path
infer_lbl <- function(p) {
  if (grepl("raj", p, ignore.case = TRUE))     return("Rajkumar (APOLLO4C)")
  if (grepl("mertins", p, ignore.case = TRUE)) return("Mertins (CPTAC BRCA)")
  if (grepl("krug", p, ignore.case = TRUE))    return("Krug (CPTAC BRCA)")
  NA_character_
}

bio_order <- c("Rajkumar (APOLLO4C)", "Mertins (CPTAC BRCA)", "Krug (CPTAC BRCA)")

# ---------------- PARSE & COLLECT ----------------
audit <- list()

# per-regime sets: list(regime -> list(UP=list(lbl=vec), DOWN=list(lbl=vec)))
sets <- list(
  FDR_only = list(UP = list(), DOWN = list()),
  FDR_FC   = list(UP = list(), DOWN = list())
)

for (p in cand_files) {
  lbl <- infer_lbl(p)
  if (is.na(lbl)) next
  
  ok <- TRUE; df <- NULL
  tryCatch({ df <- read_tab(p) }, error = function(e) ok <<- FALSE)
  if (!ok || is.null(df) || !nrow(df)) {
    audit[[p]] <- data.frame(study = lbl, file = p, used = FALSE, reason = "read_fail")
    next
  }
  if (!looks_like_dep_table(df)) {
    audit[[p]] <- data.frame(study = lbl, file = p, used = FALSE, reason = "not_dep_like")
    next
  }
  
  s1 <- make_sets_from_de(df, fdr_max = FDR_MAX, lfc_min = NULL)     # FDR-only
  s2 <- make_sets_from_de(df, fdr_max = FDR_MAX, lfc_min = LFC_MIN)  # FDR+FC
  
  # accumulate (union within cohort in case multiple tables exist)
  if (length(s1$up))   sets$FDR_only$UP[[lbl]]   <- unique(c(sets$FDR_only$UP[[lbl]],   s1$up))
  if (length(s1$down)) sets$FDR_only$DOWN[[lbl]] <- unique(c(sets$FDR_only$DOWN[[lbl]], s1$down))
  if (length(s2$up))   sets$FDR_FC$UP[[lbl]]     <- unique(c(sets$FDR_FC$UP[[lbl]],     s2$up))
  if (length(s2$down)) sets$FDR_FC$DOWN[[lbl]]   <- unique(c(sets$FDR_FC$DOWN[[lbl]],   s2$down))
  
  audit[[p]] <- data.frame(
    study   = lbl,
    file    = p,
    used    = TRUE,
    n_total = nrow(df),
    n_sig_fdr_only = s1$n_keep,
    n_sig_fdr_fc   = s2$n_keep,
    n_up_fdr_only  = length(s1$up),
    n_dn_fdr_only  = length(s1$down),
    n_up_fdr_fc    = length(s2$up),
    n_dn_fdr_fc    = length(s2$down),
    stringsAsFactors = FALSE
  )
}

audit_df <- data.table::rbindlist(audit, use.names = TRUE, fill = TRUE)
write.csv(audit_df, file.path(OUT_DIR, "PROT_overlap_audit_by_file.csv"), row.names = FALSE)

# enforce cohort ordering + drop missing cohorts
order_sets <- function(x) {
  x <- x[vapply(x, length, 1L) > 0]
  x <- x[intersect(bio_order, names(x))]
  x
}

for (rg in names(sets)) {
  sets[[rg]]$UP   <- order_sets(sets[[rg]]$UP)
  sets[[rg]]$DOWN <- order_sets(sets[[rg]]$DOWN)
}

# sanity: at least one regime has sets
if (length(sets$FDR_only$UP) + length(sets$FDR_only$DOWN) +
    length(sets$FDR_FC$UP)   + length(sets$FDR_FC$DOWN) == 0) {
  stop("No usable proteomics DE tables found. Check file patterns and column names.")
}

# ---------------- PLOTTING: UpSet ----------------
use_complex <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize", quietly = TRUE)

make_upset_complex <- function(sets_list, title_base, file_stub) {
  if (!length(sets_list)) return(invisible(NULL))
  suppressPackageStartupMessages({ library(ComplexHeatmap); library(grid); library(ragg) })
  
  cm <- ComplexHeatmap::make_comb_mat(sets_list)
  if (length(ComplexHeatmap::comb_size(cm)) == 0) return(invisible(NULL))
  
  top_idx <- order(ComplexHeatmap::comb_size(cm), decreasing = TRUE)
  top_idx <- top_idx[seq_len(min(30, length(top_idx)))]
  
  ht <- ComplexHeatmap::UpSet(
    cm[top_idx],
    pt_size = unit(3, "mm"),
    lwd = 1.2,
    top_annotation = ComplexHeatmap::upset_top_annotation(
      cm[top_idx], add_numbers = TRUE, gp = grid::gpar(fill = "gray25")
    ),
    set_order = seq_along(ComplexHeatmap::set_name(cm))
  )
  
  ragg::agg_png(file.path(OUT_DIR, paste0(file_stub, ".png")),
                width = 2600, height = 1400, units = "px", res = 200, background = "white")
  grid::grid.newpage()
  grid::grid.text(title_base, y = unit(1, "npc") - unit(8, "pt"))
  ComplexHeatmap::draw(ht, padding = unit(c(10,10,10,10), "pt"))
  dev.off()
  
  pdf(file.path(OUT_DIR, paste0(file_stub, ".pdf")), width = 17, height = 10)
  grid::grid.newpage()
  grid::grid.text(title_base, y = unit(1, "npc") - unit(8, "pt"))
  ComplexHeatmap::draw(ht, padding = unit(c(10,10,10,10), "pt"))
  dev.off()
}

make_upset_fallback <- function(sets_list, title_base, file_stub) {
  if (!length(sets_list)) return(invisible(NULL))
  studies <- names(sets_list)
  
  cc <- data.frame(combo = character(), size = integer(), order = integer(), stringsAsFactors = FALSE)
  for (s in studies) {
    cc <- rbind(cc, data.frame(combo = s, size = length(sets_list[[s]]), order = 1))
  }
  for (r in 2:length(studies)) {
    combs <- combn(studies, r, simplify = FALSE)
    for (cmb in combs) {
      inter <- Reduce(intersect, sets_list[cmb])
      cc <- rbind(cc, data.frame(combo = paste(cmb, collapse = " ∩ "),
                                 size = length(inter), order = r))
    }
  }
  cc <- cc[order(-cc$size, cc$order), ]
  topK <- head(cc, 15)
  
  png(file.path(OUT_DIR, paste0(file_stub, ".png")), width = 1700, height = 750, res = 170, bg = "white")
  op <- par(mar = c(12,5,4,1))
  barplot(topK$size, names.arg = topK$combo, las = 2,
          ylab = "Intersection size", main = title_base)
  par(op); dev.off()
  
  pdf(file.path(OUT_DIR, paste0(file_stub, ".pdf")), width = 12, height = 6)
  op <- par(mar = c(12,5,4,1))
  barplot(topK$size, names.arg = topK$combo, las = 2,
          ylab = "Intersection size", main = title_base)
  par(op); dev.off()
}

make_upset <- function(sets_list, title_base, file_stub) {
  if (use_complex) make_upset_complex(sets_list, title_base, file_stub) else make_upset_fallback(sets_list, title_base, file_stub)
}

# ---------------- PLOTTING: 3-way Venn ----------------
make_venn <- function(sets_list, title_base, file_stub) {
  if (length(sets_list) != 3) return(invisible(NULL))
  if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
  suppressPackageStartupMessages({ library(VennDiagram); library(grid) })
  
  # VennDiagram wants a named list of vectors
  lst <- sets_list
  # avoid NA/empty
  if (any(vapply(lst, length, 1L) == 0)) return(invisible(NULL))
  
  png(file.path(OUT_DIR, paste0(file_stub, ".png")), width = 1400, height = 1100, res = 160, bg = "white")
  grid::grid.newpage()
  VennDiagram::draw.triple.venn(
    area1 = length(lst[[1]]),
    area2 = length(lst[[2]]),
    area3 = length(lst[[3]]),
    n12   = length(intersect(lst[[1]], lst[[2]])),
    n23   = length(intersect(lst[[2]], lst[[3]])),
    n13   = length(intersect(lst[[1]], lst[[3]])),
    n123  = length(Reduce(intersect, lst)),
    category = names(lst),
    fill = c("gray85","gray70","gray55"),
    lty  = "solid",
    cex = 1.1,
    cat.cex = 1.2,
    cat.pos = c(-20, 20, 0),
    cat.dist = c(0.06, 0.06, 0.06),
    main = title_base
  )
  dev.off()
  
  pdf(file.path(OUT_DIR, paste0(file_stub, ".pdf")), width = 9.5, height = 8.0)
  grid::grid.newpage()
  VennDiagram::draw.triple.venn(
    area1 = length(lst[[1]]),
    area2 = length(lst[[2]]),
    area3 = length(lst[[3]]),
    n12   = length(intersect(lst[[1]], lst[[2]])),
    n23   = length(intersect(lst[[2]], lst[[3]])),
    n13   = length(intersect(lst[[1]], lst[[3]])),
    n123  = length(Reduce(intersect, lst)),
    category = names(lst),
    fill = c("gray85","gray70","gray55"),
    lty  = "solid",
    cex = 1.1,
    cat.cex = 1.2,
    cat.pos = c(-20, 20, 0),
    cat.dist = c(0.06, 0.06, 0.06),
    main = title_base
  )
  dev.off()
}

# ---------------- RUN: BOTH REGIMES + UP/DOWN ----------------
for (rg in names(sets)) {
  rg_title <- if (rg == "FDR_only") {
    sprintf("Proteome overlap — FDR≤%.2f", FDR_MAX)
  } else {
    sprintf("Proteome overlap — FDR≤%.2f & |log2FC|≥%.2f", FDR_MAX, LFC_MIN)
  }
  
  # UP
  make_upset(
    sets[[rg]]$UP,
    paste0(rg_title, " — Up in HER2+"),
    paste0("PROT_", rg, "_UpSet_UP_HER2pos")
  )
  make_venn(
    sets[[rg]]$UP,
    paste0(rg_title, " — Up in HER2+"),
    paste0("PROT_", rg, "_Venn_UP_HER2pos")
  )
  
  # DOWN
  make_upset(
    sets[[rg]]$DOWN,
    paste0(rg_title, " — Down in HER2+"),
    paste0("PROT_", rg, "_UpSet_DOWN_HER2pos")
  )
  make_venn(
    sets[[rg]]$DOWN,
    paste0(rg_title, " — Down in HER2+"),
    paste0("PROT_", rg, "_Venn_DOWN_HER2pos")
  )
}

cat("Proteome overlap UpSet/Venn outputs → ", OUT_DIR, "\n", sep = "")
