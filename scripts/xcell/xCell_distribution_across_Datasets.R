# =========================
# xCell_distribution_across_Datasets.R  
# =========================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -------- paths --------
base_dir <- "~/Desktop/Thesis/AIM 2 HER2 DATA"

files <- list(
  Wolf           = file.path(base_dir, "wolf et al",            "aim2_outputs/xcell_admixture", "Wolf_xCell_admixture_scores.tsv"),
  Robinson       = file.path(base_dir, "robinson et al",        "aim2_outputs/xcell_admixture", "Robinson_xCell_admixture_scores.tsv"),
  RajKumar       = file.path(base_dir, "raj kumar et al",       "aim2_outputs/xcell_admixture", "RK_xCell_admixture_scores.tsv"),
  Brueffer       = file.path(base_dir, "brueffer et al",        "aim2_outputs/xcell_admixture", "Brueffer_xCell_admixture_scores.tsv"),
  BarronGallardo = file.path(base_dir, "barron-gallardo et al", "aim2_outputs/xcell_admixture", "GSE163882_xCell_admixture_scores.tsv"),
  Krug           = file.path(base_dir, "krug et al",            "aim2_outputs/xcell_admixture", "Krug_xCell_admixture_scores.tsv")
)

required_cols <- c("HER2", "ImmuneScore", "StromaScore",
                   "MicroenvironmentScore", "Fibroblasts")

all_long <- list()

# -------- read and reshape each dataset --------
for (nm in names(files)) {
  path <- files[[nm]]
  message("Reading ", nm, " from: ", path)
  if (!file.exists(path)) stop("File not found: ", path)
  
  dt <- fread(path, check.names = FALSE)
  
  missing <- setdiff(required_cols, names(dt))
  if (length(missing) > 0)
    stop("Missing columns in ", nm, ": ", paste(missing, collapse = ", "))
  
  # epithelial column may be 'Epithelial.cells' or 'Epithelial cells'
  epi_col <- if ("Epithelial.cells" %in% names(dt)) {
    "Epithelial.cells"
  } else if ("Epithelial cells" %in% names(dt)) {
    "Epithelial cells"
  } else {
    stop("No epithelial column found in ", nm)
  }
  
  # reshape to long format (use data.table + then coerce)
  long <- melt(
    dt,
    id.vars      = "HER2",
    measure.vars = c("ImmuneScore", "StromaScore",
                     "MicroenvironmentScore", "Fibroblasts", epi_col),
    variable.name = "score_type",
    value.name    = "score"
  )
  long <- as.data.table(long)
  
  # unify epithelial naming
  long[, score_type := as.character(score_type)]
  long[score_type == epi_col, score_type := "Epithelial cells"]
  
  long[, dataset := nm]
  
  # consistent HER2 factor
  long[, HER2 := factor(HER2, levels = c("HER2neg", "HER2pos"))]
  
  all_long[[nm]] <- long
}

all_long <- rbindlist(all_long, use.names = TRUE, fill = TRUE)

# order facets and datasets to match slide
all_long[, score_type := factor(
  score_type,
  levels = c("ImmuneScore", "StromaScore",
             "MicroenvironmentScore", "Fibroblasts",
             "Epithelial cells")
)]

all_long[, dataset := factor(
  dataset,
  levels = c("BarronGallardo", "Brueffer", "Krug",
             "RajKumar", "Robinson", "Wolf")
)]

# -------- plot --------
COLS <- c(HER2neg = "#F8766D", HER2pos = "#00BFC4")

p <- ggplot(all_long,
            aes(x = dataset, y = score, fill = HER2)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 0.7, alpha = 0.4, color = "grey30") +
  facet_wrap(~ score_type, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = COLS, name = "HER2 group") +
  labs(
    title = "xCell Admixture Score Distributions Across Datasets",
    x = "Dataset",
    y = "xCell score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(size = 12, face = "bold"),
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold")
  )

# -------- save --------
outdir <- file.path(base_dir, "xCell_overview")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

outfile <- file.path(outdir,
                     "xCell_Admixture_Score_Distributions_Across_Datasets.png")

ggsave(outfile, p, width = 12, height = 6, dpi = 300)

message("Saved plot to: ", outfile)
