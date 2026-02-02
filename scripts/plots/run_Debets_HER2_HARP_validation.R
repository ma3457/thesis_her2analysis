# ================= run_Debets_HER2_HARP_validation.R — FINAL =================
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# -------------------------
# 1. Resolve paths
# -------------------------
here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) dirname(normalizePath(f))
  else if (!is.null(sys.frames()[[1]]$ofile)) dirname(normalizePath(sys.frames()[[1]]$ofile))
  else normalizePath(".")
})

base_debets <- here
base_root   <- dirname(base_debets)

message("Debets folder: ", base_debets)
message("Root folder  : ", base_root)

# -------------------------
# 2. PRIORITIZED FEATURE FILE — 
# -------------------------
prior_fp <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/wolf et al/aim2_outputs/Wolf_HARP_RNA/Wolf_core_HER2_activation_features_RNA_directional.csv"

if (!file.exists(prior_fp)) {
  stop("❌ PRIORITIZED FEATURE FILE NOT FOUND:\n", prior_fp)
}

# -------------------------
# 3. Debets DE output 
# -------------------------
debets_de_fp <- file.path(
  base_debets, "aim2_outputs", "Debets2023_protein", "DE_full.tsv"
)

if (!file.exists(debets_de_fp)) {
  stop("❌ Missing Debets DE_full.tsv at:\n  ", debets_de_fp)
}

out_dir <- file.path(base_debets, "aim2_outputs", "Debets_HER2_HARP_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 4. Load HER2/HARP feature list
# -------------------------
her2_feats <- fread(prior_fp)
coln <- names(her2_feats)[1]   # first column is gene name

setnames(her2_feats, coln, "gene_symbol")
her2_feats[, gene_symbol := toupper(trimws(gene_symbol))]
her2_feats <- unique(her2_feats)

message("Loaded ", nrow(her2_feats), " HER2 activation features.")

# -------------------------
# 5. Load Debets DE results
# -------------------------
deb <- fread(debets_de_fp)

# Identify columns
if (!"gene_symbol" %in% names(deb)) {
  setnames(deb, names(deb)[1], "gene_symbol")
}

lfc_col <- names(deb)[grepl("log2FC", names(deb), ignore.case = TRUE)][1]
setnames(deb, lfc_col, "log2FC_pCR_vs_No_pCR", skip_absent = TRUE)

deb[, gene_symbol := toupper(trimws(gene_symbol))]

# -------------------------
# 6. Merge + classify significance
# -------------------------
merged <- merge(
  her2_feats,
  deb,
  by = "gene_symbol",
  all.x = TRUE
)

merged[, detected_in_Debets := !is.na(log2FC_pCR_vs_No_pCR)]

# Significance classification
merged[, sig_class := fifelse(
  !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(log2FC_pCR_vs_No_pCR) >= 1,
  "FDR<0.05 & |log2FC|>=1",
  fifelse(
    !is.na(adj.P.Val) & adj.P.Val < 0.05,
    "FDR<0.05 only",
    fifelse(
      !is.na(P.Value) & P.Value < 0.05,
      "P<0.05 only",
      "NS"
    )
  )
)]

# -------------------------
# 7. Write outputs
# -------------------------
val_fp <- file.path(out_dir, "Debets_HER2_HARP_validation_table.tsv")
fwrite(merged, val_fp, sep="\t")

sig_fp <- file.path(out_dir, "Debets_HER2_HARP_significant_candidates.tsv")
fwrite(merged[sig_class != "NS" & detected_in_Debets == TRUE][order(adj.P.Val)],
       sig_fp, sep="\t")

# -------------------------
# 8. Summary
# -------------------------
total_prior   <- nrow(her2_feats)
detected      <- merged[detected_in_Debets == TRUE, .N]
any_signal    <- merged[sig_class != "NS", .N]
strong_signal <- merged[sig_class == "FDR<0.05 & |log2FC|>=1", .N]

summary_lines <- c(
  sprintf("Total prioritized HER2/HARP features: %d", total_prior),
  sprintf("Detected in Debets proteome: %d (%.1f%%)",
          detected, 100 * detected / total_prior),
  sprintf("Any signal (P<0.05 or FDR<0.05): %d", any_signal),
  sprintf("Strong signal (FDR<0.05 & |log2FC|>=1): %d", strong_signal),
  "",
  paste("Validation table:", val_fp),
  paste("Significant genes:", sig_fp)
)

writeLines(summary_lines, file.path(out_dir, "Debets_HER2_HARP_validation_summary.txt"))

message("=== DONE ===")
cat(paste(summary_lines, collapse="\n"), "\n")
# ============================================================================
