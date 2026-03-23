# consensus_HER2_transcripts.R — FINAL (robust fread + spaces-safe)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ----------------------------
# CONFIG: where your .txt lists live
# ----------------------------
base_dir <- normalizePath(
  "~/Desktop/Thesis/AIM 2 HER2 DATA/aim2_upset_outputs/transcript_upset_dual/FDR_FC",
  mustWork = FALSE
)

# If you actually saved them somewhere else, set base_dir to that folder.
# You can confirm by: list.files(base_dir)

files <- list(
  Wolf_UP       = file.path(base_dir, "GSE81538_UP_FDR0.05_LFC1.txt"),
  Wolf_DOWN     = file.path(base_dir, "GSE81538_DOWN_FDR0.05_LFC1.txt"),
  Robinson_UP   = file.path(base_dir, "GSE194040_UP_FDR0.05_LFC1.txt"),
  Robinson_DOWN = file.path(base_dir, "GSE194040_DOWN_FDR0.05_LFC1.txt"),
  Brueffer_UP   = file.path(base_dir, "GSE199633_UP_FDR0.05_LFC1.txt"),
  Brueffer_DOWN = file.path(base_dir, "GSE199633_DOWN_FDR0.05_LFC1.txt")
)

# ----------------------------
# Guardrails
# ----------------------------
missing <- names(files)[!file.exists(files)]
if (length(missing)) {
  cat("Missing files:\n")
  for (nm in missing) cat(" - ", nm, " -> ", files[[nm]], "\n", sep = "")
  stop("Fix base_dir / filenames above, then rerun.")
}

# ----------------------------
# Helpers
# ----------------------------
clean_gene <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- x[x != "" & !is.na(x)]
  x <- gsub("\\.\\d+$", "", x)     # remove ".1" suffixes
  toupper(x)
}

read_gene_list <- function(path, cohort, direction) {
  path <- normalizePath(path, mustWork = TRUE)
  
  # IMPORTANT: use file= explicitly so spaces don't trigger "system command" behavior
  dt <- data.table::fread(file = path, header = FALSE, sep = "\t", fill = TRUE, data.table = TRUE)
  
  if (ncol(dt) == 0 || nrow(dt) == 0) {
    stop("Gene list read as empty: ", path)
  }
  
  # If a header slipped in, dt[1,1] might be "gene" etc — we’ll just clean and drop non-genes later
  genes <- clean_gene(dt[[1]])
  
  data.frame(
    gene_symbol = genes,
    cohort = cohort,
    direction = direction,
    stringsAsFactors = FALSE
  )
}

# ----------------------------
# Read lists (Wolf / Robinson / Brueffer)
# ----------------------------
gene_lists <- bind_rows(
  read_gene_list(files$Wolf_UP,       "Wolf",      "UP"),
  read_gene_list(files$Wolf_DOWN,     "Wolf",      "DOWN"),
  read_gene_list(files$Robinson_UP,   "Robinson",  "UP"),
  read_gene_list(files$Robinson_DOWN, "Robinson",  "DOWN"),
  read_gene_list(files$Brueffer_UP,   "Brueffer",  "UP"),
  read_gene_list(files$Brueffer_DOWN, "Brueffer",  "DOWN")
) %>%
  distinct(gene_symbol, cohort, direction)

# ----------------------------
# Build consensus table
# ----------------------------
consensus <- gene_lists %>%
  mutate(present = TRUE) %>%
  tidyr::pivot_wider(
    names_from  = cohort,
    values_from = present,
    values_fill = list(present = FALSE)
  ) %>%
  mutate(
    n_cohorts = as.integer(Wolf) + as.integer(Robinson) + as.integer(Brueffer)
  ) %>%
  arrange(desc(n_cohorts), direction, gene_symbol)

# ----------------------------
# Write output
# ----------------------------
out_file <- file.path(base_dir, "Consensus_HER2_transcript_features.tsv")
data.table::fwrite(consensus, out_file, sep = "\t")

# ----------------------------
# Quick sanity prints
# ----------------------------
cat("Wrote: ", out_file, "\n", sep = "")
cat("Counts by n_cohorts:\n")
print(table(consensus$n_cohorts))

cat("\nCanonical check:\n")
print(consensus %>% filter(gene_symbol %in% c("ERBB2","GRB7","STARD3","MIEN1")))
