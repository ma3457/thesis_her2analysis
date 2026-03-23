library(dplyr)
library(ggplot2)

# Force join keys to be the same type
p76p25 <- p76p25 %>% mutate(SampleID = as.character(SampleID))
harps  <- harps  %>% mutate(SampleID = as.character(SampleID))
meta   <- meta   %>% mutate(SampleID = as.character(SampleID))

# Join (keep everything from p76p25 as the base)
joint <- p76p25 %>%
  left_join(harps, by = "SampleID") %>%
  left_join(meta,  by = "SampleID")

# Quick sanity checks
cat("Rows in p76p25:", nrow(p76p25), "\n")
cat("Rows in joint :", nrow(joint), "\n")
cat("Missing HARPS:", sum(is.na(joint$HARPS_score)), "\n")
cat("Missing Group:", sum(is.na(joint$Group)), "\n")

# Plot: P76 vs HARPS (color by BP HER2 group from meta_USED)
p <- ggplot(joint, aes(x = P76_meanZ, y = HARPS_score, color = Group)) +
  geom_point(alpha = 0.75) +
  theme_classic() +
  labs(
    x = "P76 mean Z score",
    y = "HARPS score",
    title = "Wolf: Joint pathway score space (P76 vs HARPS)"
  )

print(p)

ggsave(file.path(out_dir, "joint_scores_scatter_P76meanZ_vs_HARPS.png"),
       plot = p, width = 7.5, height = 6.0, dpi = 180)

# Optional second plot: P25 vs HARPS
p2 <- ggplot(joint, aes(x = P25_meanZ, y = HARPS_score, color = Group)) +
  geom_point(alpha = 0.75) +
  theme_classic() +
  labs(
    x = "P25 mean Z score",
    y = "HARPS score",
    title = "Wolf: Joint pathway score space (P25 vs HARPS)"
  )

print(p2)

ggsave(file.path(out_dir, "joint_scores_scatter_P25meanZ_vs_HARPS.png"),
       plot = p2, width = 7.5, height = 6.0, dpi = 180)

# Save the merged table
fwrite(as.data.table(joint), file.path(out_dir, "joint_scores_merged.tsv"), sep = "\t")

library(ggplot2)

# 1) Money plot: P76 vs HARPS, colored by HER2 group
p <- ggplot(joint, aes(x = P76_meanZ, y = HARPS_score, color = Group)) +
  geom_point(alpha = 0.75) +
  theme_classic() +
  labs(
    x = "P76 mean Z score",
    y = "HARPS score",
    title = "Wolf: Joint pathway score space (P76 vs HARPS)"
  )

print(p)

ggsave(
  filename = file.path(out_dir, "FIG_joint_P76_vs_HARPS_Wolf.png"),
  plot = p, width = 7.5, height = 6.0, dpi = 220
)

# 2) One number to report: Spearman correlation
rho <- cor(joint$P76_meanZ, joint$HARPS_score, method = "spearman", use = "pairwise.complete.obs")
rho
writeLines(
  sprintf("Wolf Spearman rho(P76_meanZ, HARPS_score) = %.3f", rho),
  con = file.path(out_dir, "FIG_joint_P76_vs_HARPS_Wolf_Spearman.txt")
)
