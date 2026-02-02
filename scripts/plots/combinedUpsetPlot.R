# Aim2UpSet_Combined_roots.R — Combined UpSet

suppressPackageStartupMessages({
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  if (!requireNamespace("ragg", quietly = TRUE))        install.packages("ragg")
  library(data.table)
})

# -------------------------- CONFIG -------------------------------------------
FDR_MAX <- 0.05
LFC_MIN <- 1.00   # FC2

base_dir <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA"
out_root <- file.path(base_dir, "aim2_upset_outputs")
out_dir  <- file.path(out_root, "combined_upset")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------- ROOTS -------------------------------

# Transcript HER2+ vs HER2−
tx_roots <- c(
  "GSE194040_ISPY2_mRNA_BPsubtype"   =
    file.path(base_dir, "wolf et al/aim2_outputs/GSE194040_ISPY2_mRNA_BPsubtype"),
  "GSE199633_Robinson2025"           =
    file.path(base_dir, "robinson et al/aim2_outputs/GSE199633_Robinson2025"),
  "GSE293591_Kushnarev2025_ERBB2cut" =
    file.path(base_dir, "kushnarev et al/aim2_outputs/GSE293591_Kushnarev2025_ERBB2cut"),
  "GSE212143_Greer2022_celltype"     =
    file.path(base_dir, "greer et al/aim2_outputs/GSE212143_Greer2022_celltype"),
  "GSE81538"                         =
    file.path(base_dir, "brueffer et al/aim2_outputs/GSE81538"),
  "GSE163882"                        =
    file.path(base_dir, "barron-gallardo et al/aim2_outputs/GSE163882"),
  "Krug_CPTAC_BRCA_RNA_A"           =
    file.path(base_dir, "krug et al/aim2_outputs/Krug_CPTAC_BRCA_RNA_RNA_A_ClinicalOnly"),
  "RajKumar_RNA"                     =
    file.path(base_dir, "raj kumar et al/aim2_outputs/RajKumar_RNA")
)

# Proteome / RPPA roots (we'll scan inside these)
prot_roots <- c(
  Wolf_ISPY2_RPPA  = file.path(base_dir, "wolf et al/aim2_outputs"),
  Debets2023       = file.path(base_dir, "Debets et al/aim2_outputs"),
  RajKumar         = file.path(base_dir, "raj kumar et al/aim2_outputs"),
  Krug_CPTAC_BRCA  = file.path(base_dir, "krug et al/aim2_outputs")
)

# pretty labels for combined UpSet
label_map <- c(
  "GSE194040_ISPY2_mRNA_BPsubtype"   = "Wolf (Transcript)",
  "GSE199633_Robinson2025"           = "Robinson (Transcript)",
  "GSE293591_Kushnarev2025_ERBB2cut" = "Kushnarev (Transcript)",
  "GSE212143_Greer2022_celltype"     = "Greer (Transcript)",
  "GSE81538"                         = "Brueffer (Transcript)",
  "GSE163882"                        = "Barron–Gallardo (Tx)",
  "Krug_CPTAC_BRCA_RNA_A"            = "Krug CPTAC (Transcript)",
  "RajKumar_RNA"                     = "RajKumar (Transcript)",
  "Wolf_ISPY2_RPPA"                  = "Wolf (RPPA)",
  "Debets2023"                       = "Debets (Proteome)",
  "RajKumar"                         = "RajKumar (Proteome)",
  "Krug_CPTAC_BRCA"                  = "Krug CPTAC (Proteome)"
)

# biological row order
combined_order <- c(
  "Wolf (Transcript)",
  "Robinson (Transcript)",
  "Kushnarev (Transcript)",
  "Greer (Transcript)",
  "Brueffer (Transcript)",
  "Barron–Gallardo (Tx)",
  "Krug CPTAC (Transcript)",
  "RajKumar (Transcript)",
  "Wolf (RPPA)",
  "Krug CPTAC (Proteome)",
  "Debets (Proteome)",
  "RajKumar (Proteome)"
)

# -------------------------- SHARED HELPERS -----------------------------------

read_tab <- function(p) {
  if (grepl("\\.csv$", p, TRUE)) read.csv(p, check.names = FALSE)
  else read.table(p, header = TRUE, sep = "\t", quote = "",
                  comment.char = "", check.names = FALSE)
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

clean_syms_tx <- function(v) {
  v <- as.character(v)
  v <- gsub("^\\s+|\\s+$", "", v)
  v <- gsub("\\.\\d+$", "", v)
  toupper(v)
}

clean_syms_prot <- function(v) {
  v <- as.character(v)
  v <- gsub("^\\s+|\\s+$", "", v)
  v <- gsub("\\.\\d+$", "", v)                      # version
  v <- gsub("\\.total$", "", v, ignore.case=TRUE)   # RPPA suffix
  v <- gsub("\\.pY\\d+$", "", v, ignore.case=TRUE)  # phospho tag
  v <- gsub("_TOTAL$", "", v, ignore.case=TRUE)
  v <- gsub("[- ]", "_", v)
  toupper(v)
}

# simple: any file with HER2pos_vs_HER2neg in name
pick_de_file_tx <- function(pathlike) {
  if (file.exists(pathlike) && !dir.exists(pathlike)) return(pathlike)
  if (!dir.exists(pathlike)) return(NA_character_)
  cand <- list.files(
    pathlike,
    pattern = "HER2pos_vs_HER2neg.*\\.(tsv|csv)$",
    full.names = TRUE, ignore.case = TRUE, recursive = TRUE
  )
  if (!length(cand)) return(NA_character_)
  # prefer ALL_* then FILT_* if both exist
  cand_all  <- grep("^ALL_.*HER2pos_vs_HER2neg", basename(cand), value = TRUE)
  if (length(cand_all)) {
    cand <- cand[match(cand_all[1], basename(cand))]
  } else {
    cand <- cand[1]
  }
  cand
}

# proteome-ish detection from your proteome script
is_proteomeish_path <- function(p) {
  grepl("RPPA|PROT|PROTEOM|protein", p, ignore.case = TRUE)
}
looks_like_dep_table <- function(df) {
  cn <- tolower(names(df))
  has_fc  <- any(c("logfc","log2foldchange","log2fc","log2fc_her2pos_vs_her2neg") %in% cn)
  has_fdr <- any(c("adj.p.val","padj","fdr","qval","adj_pval","fdr.bh") %in% cn)
  has_id  <- any(c("gene_symbol","symbol","gene","genesymbol","feature_id","protein") %in% cn)
  has_fc && has_fdr && has_id
}

# -------------------------- TRANSCRIPTS --------------------------------------

tx_audit <- list(); tx_sets_up <- list(); tx_sets_down <- list()

for (id in names(tx_roots)) {
  root <- tx_roots[[id]]
  p <- pick_de_file_tx(root)
  if (is.na(p)) {
    tx_audit[[id]] <- data.frame(study=id, file=NA, n_total=NA,
                                 n_sig=NA, n_up=NA, n_down=NA,
                                 note="no_HER2_DE_file")
    next
  }
  dt <- read_tab(p)
  g <- pick_col(names(dt),
                c("gene_symbol","SYMBOL","symbol","Gene","GeneSymbol","gene","feature_id"))
  q <- pick_col(names(dt),
                c("adj.P.Val","padj","FDR","qval","adj_pval","FDR.BH"))
  l <- pick_col(names(dt),
                c("logFC","log2FoldChange","log2FC","log2FC_HER2pos_vs_HER2neg"))
  if (any(is.na(c(g,q,l)))) {
    tx_audit[[id]] <- data.frame(study=id, file=p, n_total=nrow(dt),
                                 n_sig=NA, n_up=NA, n_down=NA,
                                 note="missing_cols")
    next
  }
  fdr <- suppressWarnings(as.numeric(dt[[q]]))
  lfc <- suppressWarnings(as.numeric(dt[[l]]))
  sym <- clean_syms_tx(dt[[g]])
  ok  <- which(!is.na(fdr) & !is.na(lfc) & !is.na(sym) &
                 fdr <= FDR_MAX & abs(lfc) >= LFC_MIN)
  up_syms   <- unique(sym[ok][lfc[ok] > 0])
  down_syms <- unique(sym[ok][lfc[ok] < 0])
  
  tx_sets_up[[id]]   <- up_syms
  tx_sets_down[[id]] <- down_syms
  
  tx_audit[[id]] <- data.frame(
    study=id, file=p, n_total=nrow(dt),
    n_sig=length(ok), n_up=length(up_syms),
    n_down=length(down_syms), note=""
  )
}

# -------------------------- PROTEOME / RPPA ----------------------------------

prot_audit <- list(); prot_sets_up <- list(); prot_sets_down <- list()

for (rid in names(prot_roots)) {
  base <- prot_roots[[rid]]
  if (!dir.exists(base)) {
    prot_audit[[rid]] <- data.frame(study=rid, file=NA, used=FALSE, reason="root_missing")
    next
  }
  # candidate files: proteome-ish paths + HER2pos_vs_HER2neg in name
  cand <- list.files(base, pattern = "\\.(tsv|csv)$",
                     full.names = TRUE, recursive = TRUE)
  cand <- cand[is_proteomeish_path(cand) | is_proteomeish_path(dirname(cand))]
  cand <- cand[grepl("HER2pos_vs_HER2neg", cand, ignore.case = TRUE)]
  cand <- unique(cand)
  if (!length(cand)) {
    prot_audit[[rid]] <- data.frame(study=rid, file=NA, used=FALSE, reason="no_HER2_prot_files")
    next
  }
  used_any <- FALSE
  for (p in cand) {
    ok <- TRUE; dt <- NULL
    tryCatch({ dt <- read_tab(p) }, error = function(e) ok <<- FALSE)
    if (!ok || is.null(dt) || !nrow(dt) || !looks_like_dep_table(dt)) {
      prot_audit[[paste0(rid,"|",basename(p))]] <-
        data.frame(study=rid, file=p, used=FALSE, reason="not_dep_like")
      next
    }
    g <- pick_col(names(dt),
                  c("gene_symbol","SYMBOL","symbol","Gene","GeneSymbol",
                    "gene","feature_id","Protein"))
    q <- pick_col(names(dt),
                  c("adj.P.Val","padj","FDR","qval","adj_pval","FDR.BH"))
    l <- pick_col(names(dt),
                  c("logFC","log2FoldChange","log2FC","log2FC_HER2pos_vs_HER2neg"))
    if (any(is.na(c(g,q,l)))) {
      prot_audit[[paste0(rid,"|",basename(p))]] <-
        data.frame(study=rid, file=p, used=FALSE, reason="missing_cols")
      next
    }
    fdr <- suppressWarnings(as.numeric(dt[[q]]))
    lfc <- suppressWarnings(as.numeric(dt[[l]]))
    syms <- clean_syms_prot(dt[[g]])
    
    keep <- which(!is.na(fdr) & !is.na(lfc) & !is.na(syms) &
                    fdr <= FDR_MAX & abs(lfc) >= LFC_MIN)
    up_syms   <- unique(syms[keep][lfc[keep] > 0])
    down_syms <- unique(syms[keep][lfc[keep] < 0])
    
    if (length(up_syms)) {
      prot_sets_up[[rid]] <- unique(c(prot_sets_up[[rid]], up_syms))
    }
    if (length(down_syms)) {
      prot_sets_down[[rid]] <- unique(c(prot_sets_down[[rid]], down_syms))
    }
    
    prot_audit[[paste0(rid,"|",basename(p))]] <- data.frame(
      study=rid, file=p, used=TRUE,
      n_total=nrow(dt), n_sig=length(keep),
      n_up=length(up_syms), n_down=length(down_syms),
      stringsAsFactors = FALSE
    )
    used_any <- TRUE
  }
  if (!used_any && !length(cand)) {
    prot_audit[[rid]] <- data.frame(study=rid, file=NA,
                                    used=FALSE, reason="no_valid_dep_tables")
  }
}

# -------------------------- MERGE + LABEL ------------------------------------

# drop empties
tx_sets_up   <- tx_sets_up[vapply(tx_sets_up, length, 1L) > 0]
tx_sets_down <- tx_sets_down[vapply(tx_sets_down, length, 1L) > 0]
prot_sets_up   <- prot_sets_up[vapply(prot_sets_up, length, 1L) > 0]
prot_sets_down <- prot_sets_down[vapply(prot_sets_down, length, 1L) > 0]

sets_up   <- c(tx_sets_up,   prot_sets_up)
sets_down <- c(tx_sets_down, prot_sets_down)

stopifnot(length(sets_up) + length(sets_down) > 0)

# relabel to pretty names
if (length(sets_up)) {
  new_names <- unname(label_map[names(sets_up)])
  names(sets_up) <- ifelse(is.na(new_names), names(sets_up), new_names)
}
if (length(sets_down)) {
  new_names <- unname(label_map[names(sets_down)])
  names(sets_down) <- ifelse(is.na(new_names), names(sets_down), new_names)
}

# enforce biological order
if (length(sets_up)) {
  present <- intersect(combined_order, names(sets_up))
  if (length(present)) sets_up <- sets_up[present]
}
if (length(sets_down)) {
  present <- intersect(combined_order, names(sets_down))
  if (length(present)) sets_down <- sets_down[present]
}

# -------------------------- AUDIT WRITE --------------------------------------

tx_audit_df   <- if (length(tx_audit))   rbindlist(tx_audit,   use.names = TRUE, fill = TRUE) else data.frame()
prot_audit_df <- if (length(prot_audit)) rbindlist(prot_audit, use.names = TRUE, fill = TRUE) else data.frame()
audit_df <- rbind(tx_audit_df, prot_audit_df, fill = TRUE)
if (nrow(audit_df)) {
  write.csv(audit_df,
            file.path(out_dir, "COMBINED_DE_audit_counts.csv"),
            row.names = FALSE)
}

# -------------------------- PLOT ---------------------------------------------

use_complex <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize", quietly = TRUE)

if (use_complex) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap); library(circlize); library(ragg); library(grid)
  })
  
  make_upset <- function(sets_list, title_base, file_stub) {
    if (!length(sets_list)) return(invisible(NULL))
    cm <- ComplexHeatmap::make_comb_mat(sets_list)
    if (length(ComplexHeatmap::comb_size(cm)) == 0) return(invisible(NULL))
    
    # keep the set order we already imposed with combined_order
    set_order <- seq_along(ComplexHeatmap::set_name(cm))
    top_idx   <- order(ComplexHeatmap::comb_size(cm), decreasing = TRUE)
    top_idx   <- top_idx[seq_len(min(30, length(top_idx)))]
    
    ht <- ComplexHeatmap::UpSet(
      cm[top_idx],
      pt_size = unit(3, "mm"), lwd = 1.2,
      top_annotation = ComplexHeatmap::upset_top_annotation(
        cm[top_idx], add_numbers = TRUE, gp = gpar(fill = "gray25")),
      set_order = set_order
    )
    
    ragg::agg_png(file.path(out_dir, paste0(file_stub, ".png")),
                  width = 2600, height = 1400, res = 200, background = "white")
    grid::grid.newpage()
    grid::grid.text(
      sprintf("%s — FDR≤%.2f, |log2FC|≥%.2f", title_base, FDR_MAX, LFC_MIN),
      y = unit(1, "npc") - unit(8, "pt")
    )
    ComplexHeatmap::draw(ht, padding = unit(c(10,10,10,10), "pt"))
    dev.off()
    
    pdf(file.path(out_dir, paste0(file_stub, ".pdf")),
        width = 17, height = 10, onefile = FALSE)
    grid::grid.newpage()
    grid::grid.text(
      sprintf("%s — FDR≤%.2f, |log2FC|≥%.2f", title_base, FDR_MAX, LFC_MIN),
      y = unit(1, "npc") - unit(8, "pt")
    )
    ComplexHeatmap::draw(ht, padding = unit(c(10,10,10,10), "pt"))
    dev.off()
  }
  
  if (length(sets_up))
    make_upset(
      sets_up,
      "Combined UpSet (Transcript + Proteome/RPPA): Up in HER2+",
      "Combined_UpSet_UP"
    )
  if (length(sets_down))
    make_upset(
      sets_down,
      "Combined UpSet (Transcript + Proteome/RPPA): Down in HER2+",
      "Combined_UpSet_DOWN"
    )
  
} else {
  # simple barplot fallback
  make_bar <- function(sets_list, title, stub) {
    if (!length(sets_list)) return(invisible(NULL))
    studies <- names(sets_list)
    cc <- data.frame(combo=character(), size=integer(), order=integer())
    for (s in studies) cc <- rbind(cc, data.frame(combo=s, size=length(sets_list[[s]]), order=1))
    for (r in 2:3) {
      combs <- combn(studies, r, simplify = FALSE)
      for (cmb in combs) {
        inter <- Reduce(intersect, sets_list[cmb])
        cc <- rbind(cc, data.frame(combo=paste(cmb, collapse=" ∩ "), size=length(inter), order=r))
      }
    }
    cc <- cc[order(-cc$size, cc$order), ]
    topK <- head(cc, 20)
    
    png(file.path(out_dir, paste0(stub, ".png")),
        width = 1600, height = 700, res = 160, bg = "white")
    op <- par(mar = c(12,5,4,1))
    barplot(
      topK$size, names.arg = topK$combo, las = 2,
      ylab = "Intersection size (sig features)",
      main = sprintf("%s — FDR≤%.2f, |log2FC|≥%.2f", title, FDR_MAX, LFC_MIN)
    )
    par(op); dev.off()
    
    pdf(file.path(out_dir, paste0(stub, ".pdf")),
        width = 12, height = 5.5)
    op <- par(mar = c(12,5,4,1))
    barplot(
      topK$size, names.arg = topK$combo, las = 2,
      ylab = "Intersection size (sig features)",
      main = sprintf("%s — FDR≤%.2f, |log2FC|≥%.2f", title, FDR_MAX, LFC_MIN)
    )
    par(op); dev.off()
  }
  
  make_bar(sets_up,   "Combined (Tx + Prot/RPPA) — Up in HER2+",   "Combined_UpSet_UP")
  make_bar(sets_down, "Combined (Tx + Prot/RPPA) — Down in HER2+", "Combined_UpSet_DOWN")
}

cat("Wrote combined UpSet outputs to: ", out_dir, "\n", sep = "")
