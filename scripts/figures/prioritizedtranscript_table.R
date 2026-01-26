# ================== prioritizedtranscript_table.R ==================
# prioritized UP-only overlap list (e.g., 76 genes)
# + per-study logFC, P.Value, adj.P.Val (FDR) for Wolf/Robinson/Brueffer.
# ALSO: derives the stricter 25-gene set (FDR + |log2FC|>=1 in ALL THREE cohorts, UP-only)

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# -------------------------------------------------------------------
# INPUTS (EDIT ONLY IF YOUR PATHS DIFFER)
# -------------------------------------------------------------------
PRIOR_UP_FILE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_upset_outputs/transcript_upset_dual_3cohorts_UPonly/FDRonly/HER2_Prioritized_UP_Overlap_Phase1_FDR0.05.csv"

# IMPORTANT: set these to your real ALL_DE outputs (full tables)
WOLF_ALL_DE     <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/GSE194040_ISPY2_mRNA_BPsubtype/ALL_HER2pos_vs_HER2neg_DE.tsv"
ROBINSON_ALL_DE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/robinson et al/aim2_outputs/GSE199633_Robinson2025/ALL_HER2pos_vs_HER2neg_DE.tsv"
BRUEFFER_ALL_DE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/brueffer et al/aim2_outputs/GSE81538/ALL_HER2pos_vs_HER2neg_DE.tsv"

OUT_DIR <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# thresholds used to derive the stricter set (P25)
FDR_MAX <- 0.05
LFC_MIN <- 1.0

OUT_FILE_76   <- file.path(OUT_DIR, "HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv")
OUT_FILE_25   <- file.path(OUT_DIR, "HER2_UP_overlap_25_FDR0.05_LFC1_ALL3_with_stats_Wolf_Robinson_Brueffer.csv")
OUT_LIST_25   <- file.path(OUT_DIR, "HER2_UP_overlap_25_gene_list.csv")
AUDIT_OUT_76  <- file.path(OUT_DIR, "HER2_UP_overlap_76_with_stats_AUDIT.txt")
AUDIT_OUT_25  <- file.path(OUT_DIR, "HER2_UP_overlap_25_with_stats_AUDIT.txt")

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
clean_gene <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- x[x != "" & !is.na(x)]
  x <- gsub("\\.\\d+$", "", x)   # strip ".1"
  toupper(x)
}

pick_first_present <- function(nms, cands) {
  hit <- intersect(cands, nms)
  if (length(hit)) return(hit[1])
  NA_character_
}

read_prior_up <- function(path) {
  stopifnot(file.exists(path))
  dt <- fread(file = path)
  gcol <- pick_first_present(names(dt), c("gene_symbol","feature_id","SYMBOL","symbol","Gene","gene"))
  if (is.na(gcol)) gcol <- names(dt)[1]
  unique(clean_gene(dt[[gcol]]))
}

read_de_standard <- function(path, label) {
  stopifnot(file.exists(path))
  dt <- fread(file = path)
  
  gene_col <- pick_first_present(names(dt), c("feature_id","gene_symbol","SYMBOL","symbol","Gene","gene"))
  if (is.na(gene_col)) stop("Cannot find gene column in: ", path)
  
  logfc_col <- pick_first_present(names(dt), c("logFC","log2FC","log2FoldChange"))
  p_col     <- pick_first_present(names(dt), c("P.Value","PValue","p.value","pval"))
  fdr_col   <- pick_first_present(names(dt), c("adj.P.Val","adj.P.Val.","padj","FDR","qval","adj_pval"))
  
  if (is.na(logfc_col)) stop("No logFC column in: ", path)
  if (is.na(p_col))     stop("No P.Value column in: ", path)
  if (is.na(fdr_col))   stop("No adj.P.Val/FDR column in: ", path)
  
  out <- data.table(
    gene_symbol = clean_gene(dt[[gene_col]]),
    logFC       = suppressWarnings(as.numeric(dt[[logfc_col]])),
    P.Value     = suppressWarnings(as.numeric(dt[[p_col]])),
    adj.P.Val   = suppressWarnings(as.numeric(dt[[fdr_col]]))
  )
  
  # de-dupe: keep most significant row per gene (ties -> larger |logFC|)
  out[, abs_logFC := abs(logFC)]
  setorder(out, gene_symbol, adj.P.Val, P.Value, -abs_logFC)
  out <- out[!duplicated(gene_symbol)]
  out[, abs_logFC := NULL]
  
  setnames(out,
           old = c("logFC","P.Value","adj.P.Val"),
           new = paste0(label, c("_logFC", "_P.Value", "_adj.P.Val")))
  out
}

dir_from_fc <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x > 0, "UP", "DOWN"))
}

# stricter filter used for the "25" (UP-only, ALL 3 cohorts pass FDR + FC)
pass_strict <- function(logfc, fdr, fdr_max = 0.05, lfc_min = 1.0) {
  !is.na(logfc) & !is.na(fdr) & (fdr < fdr_max) & (logfc >= lfc_min)
}

# -------------------------------------------------------------------
# Run
# -------------------------------------------------------------------
stopifnot(file.exists(PRIOR_UP_FILE))
prior <- read_prior_up(PRIOR_UP_FILE)

wolf     <- read_de_standard(WOLF_ALL_DE, "Wolf")
robinson <- read_de_standard(ROBINSON_ALL_DE, "Robinson")
brueffer <- read_de_standard(BRUEFFER_ALL_DE, "Brueffer")

m <- data.table(gene_symbol = prior)
m <- merge(m, wolf,     by = "gene_symbol", all.x = TRUE)
m <- merge(m, robinson, by = "gene_symbol", all.x = TRUE)
m <- merge(m, brueffer, by = "gene_symbol", all.x = TRUE)

# direction sanity (should all be UP by construction, but check anyway)
m[, Wolf_dir     := dir_from_fc(Wolf_logFC)]
m[, Robinson_dir := dir_from_fc(Robinson_logFC)]
m[, Brueffer_dir := dir_from_fc(Brueffer_logFC)]
m[, direction_ok_all_UP := (Wolf_dir == "UP" & Robinson_dir == "UP" & Brueffer_dir == "UP")]

# sort by strongest evidence (min FDR across cohorts)
m[, min_FDR := pmin(Wolf_adj.P.Val, Robinson_adj.P.Val, Brueffer_adj.P.Val, na.rm = TRUE)]
setorder(m, min_FDR, gene_symbol)

# ------------------------
# Output 76 table (stats)
# ------------------------
fwrite(m, OUT_FILE_76)

audit76 <- c(
  "HER2 prioritized UP overlap (FDR-only) with per-study stats",
  paste0("PRIOR_UP_FILE: ", PRIOR_UP_FILE),
  paste0("WOLF_ALL_DE: ", WOLF_ALL_DE),
  paste0("ROBINSON_ALL_DE: ", ROBINSON_ALL_DE),
  paste0("BRUEFFER_ALL_DE: ", BRUEFFER_ALL_DE),
  "",
  paste0("n genes in prioritized list (expected ~76): ", nrow(m)),
  paste0("missing Wolf stats: ", sum(is.na(m$Wolf_logFC))),
  paste0("missing Robinson stats: ", sum(is.na(m$Robinson_logFC))),
  paste0("missing Brueffer stats: ", sum(is.na(m$Brueffer_logFC))),
  paste0("fail all-UP direction check: ", sum(!m$direction_ok_all_UP, na.rm = TRUE)),
  "",
  paste0("Wrote: ", OUT_FILE_76)
)
writeLines(audit76, AUDIT_OUT_76)

cat("Wrote 76 table:\n  ", OUT_FILE_76, "\n", sep = "")
cat("Wrote 76 audit:\n  ", AUDIT_OUT_76, "\n", sep = "")

# ------------------------
# Derive + output 25 set
# ------------------------
m[, pass_Wolf     := pass_strict(Wolf_logFC,     Wolf_adj.P.Val,     FDR_MAX, LFC_MIN)]
m[, pass_Robinson := pass_strict(Robinson_logFC, Robinson_adj.P.Val, FDR_MAX, LFC_MIN)]
m[, pass_Brueffer := pass_strict(Brueffer_logFC, Brueffer_adj.P.Val, FDR_MAX, LFC_MIN)]

m25 <- m[pass_Wolf & pass_Robinson & pass_Brueffer]

# keep same ordering logic (min FDR) for readability
m25[, min_FDR := pmin(Wolf_adj.P.Val, Robinson_adj.P.Val, Brueffer_adj.P.Val, na.rm = TRUE)]
setorder(m25, min_FDR, gene_symbol)

fwrite(m25, OUT_FILE_25)
fwrite(data.table(gene_symbol = m25$gene_symbol), OUT_LIST_25)

audit25 <- c(
  "HER2 prioritized UP overlap (derived stricter set) with per-study stats",
  paste0("Derived from PRIOR_UP_FILE + thresholds across ALL 3 cohorts:"),
  paste0("  FDR_MAX = ", FDR_MAX),
  paste0("  LFC_MIN = ", LFC_MIN, " (UP-only, i.e., logFC >= 1)"),
  "",
  paste0("Input PRIOR_UP_FILE: ", PRIOR_UP_FILE),
  paste0("WOLF_ALL_DE: ", WOLF_ALL_DE),
  paste0("ROBINSON_ALL_DE: ", ROBINSON_ALL_DE),
  paste0("BRUEFFER_ALL_DE: ", BRUEFFER_ALL_DE),
  "",
  paste0("n genes in 76 list: ", nrow(m)),
  paste0("n genes passing strict ALL-3 filter (expected ~25): ", nrow(m25)),
  "",
  paste0("Wrote: ", OUT_FILE_25),
  paste0("Wrote gene list: ", OUT_LIST_25)
)
writeLines(audit25, AUDIT_OUT_25)

cat("Wrote 25 table:\n  ", OUT_FILE_25, "\n", sep = "")
cat("Wrote 25 gene list:\n  ", OUT_LIST_25, "\n", sep = "")
cat("Wrote 25 audit:\n  ", AUDIT_OUT_25, "\n", sep = "")

# =====================================================================
