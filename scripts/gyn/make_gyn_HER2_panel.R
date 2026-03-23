# make_gyn_HER2_panel.R
# 4-panel figure for gyn HER2 projection results

library(tidyverse)
library(patchwork)

# ----------------------------------------
# file paths
# ----------------------------------------
ovarian_scores_file <- "/Users/maya.anand/Desktop/Thesis/GYN DATA/ovarian_HER2_scores.csv"
ucec_scores_file    <- "/Users/maya.anand/Desktop/Thesis/GYN DATA/ucec_HER2_scores.csv"

output_dir  <- "/Users/maya.anand/Desktop/Thesis/figures"
output_file <- file.path(output_dir, "Figure_gyn_HER2_projection_panel.png")

# ----------------------------------------
# read score tables
# ----------------------------------------
ovarian_df <- read_csv(ovarian_scores_file, show_col_types = FALSE)
ucec_df    <- read_csv(ucec_scores_file, show_col_types = FALSE)

# ----------------------------------------
# helper for correlation labels
# ----------------------------------------
make_cor_label <- function(df, ycol) {
  test <- cor.test(df$ERBB2, df[[ycol]], method = "spearman")
  rho  <- unname(test$estimate)
  pval <- test$p.value
  
  p_text <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", signif(pval, 2)))
  paste0("Spearman \u03C1 = ", round(rho, 3), "\n", p_text)
}

ovarian_p76_label <- make_cor_label(ovarian_df, "P76_score")
ovarian_p25_label <- make_cor_label(ovarian_df, "P25_score")
ucec_p76_label    <- make_cor_label(ucec_df, "P76_score")
ucec_p25_label    <- make_cor_label(ucec_df, "P25_score")

# ----------------------------------------
# helper to place annotation nicely
# ----------------------------------------
get_label_pos <- function(df, ycol) {
  list(
    x = quantile(df$ERBB2, 0.05, na.rm = TRUE),
    y = quantile(df[[ycol]], 0.95, na.rm = TRUE)
  )
}

pos_ov_p76 <- get_label_pos(ovarian_df, "P76_score")
pos_ov_p25 <- get_label_pos(ovarian_df, "P25_score")
pos_uc_p76 <- get_label_pos(ucec_df, "P76_score")
pos_uc_p25 <- get_label_pos(ucec_df, "P25_score")

# ----------------------------------------
# common theme
# ----------------------------------------
panel_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# ----------------------------------------
# panel A - ovarian P76
# ----------------------------------------
p1 <- ggplot(ovarian_df, aes(x = ERBB2, y = P76_score)) +
  geom_point(size = 2.3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  annotate(
    "text",
    x = pos_ov_p76$x,
    y = pos_ov_p76$y,
    label = ovarian_p76_label,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  labs(
    title = "A  Ovarian cancer: P76 HER2 signature vs ERBB2 protein",
    x = "ERBB2 protein abundance (z-score)",
    y = "HER2 activation score"
  ) +
  panel_theme

# ----------------------------------------
# panel B - ovarian P25
# ----------------------------------------
p2 <- ggplot(ovarian_df, aes(x = ERBB2, y = P25_score)) +
  geom_point(size = 2.3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  annotate(
    "text",
    x = pos_ov_p25$x,
    y = pos_ov_p25$y,
    label = ovarian_p25_label,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  labs(
    title = "B  Ovarian cancer: P25 HER2 signature vs ERBB2 protein",
    x = "ERBB2 protein abundance (z-score)",
    y = "HER2 activation score"
  ) +
  panel_theme

# ----------------------------------------
# panel C - UCEC P76
# ----------------------------------------
p3 <- ggplot(ucec_df, aes(x = ERBB2, y = P76_score)) +
  geom_point(size = 2.3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  annotate(
    "text",
    x = pos_uc_p76$x,
    y = pos_uc_p76$y,
    label = ucec_p76_label,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  labs(
    title = "C  Endometrial carcinoma (UCEC): P76 HER2 signature vs ERBB2 protein",
    x = "ERBB2 protein abundance (z-score)",
    y = "HER2 activation score"
  ) +
  panel_theme

# ----------------------------------------
# panel D - UCEC P25
# ----------------------------------------
p4 <- ggplot(ucec_df, aes(x = ERBB2, y = P25_score)) +
  geom_point(size = 2.3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  annotate(
    "text",
    x = pos_uc_p25$x,
    y = pos_uc_p25$y,
    label = ucec_p25_label,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  labs(
    title = "D  Endometrial carcinoma (UCEC): P25 HER2 signature vs ERBB2 protein",
    x = "ERBB2 protein abundance (z-score)",
    y = "HER2 activation score"
  ) +
  panel_theme

# ----------------------------------------
# combine into one panel
# ----------------------------------------
final_panel <- (p1 | p2) / (p3 | p4)

# ----------------------------------------
# save figure
# ----------------------------------------
ggsave(
  filename = output_file,
  plot = final_panel,
  width = 16,
  height = 12,
  dpi = 300
)

cat("\nSaved final panel to:\n", output_file, "\n")
