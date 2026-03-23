# ================= Robinson P76/P25 Scoring =================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# --------- Paths ---------
base_dir <- "~/Desktop/Thesis/AIM 2 HER2 DATA/robinson et al"
res_dir  <- file.path(base_dir, "aim2_outputs", "GSE199633_Robinson2025")

expr_file <- file.path(base_dir, "GSE199633_expr_for_pipeline.tsv")
meta_file <- file.path(res_dir, "GSE199633_Robinson2025_RNA_meta_USED.tsv")

prior_dir <- "~/Desktop/Thesis/AIM 2 HER2 DATA/aim2_prioritized"

P76_file <- file.path(prior_dir,
                      "HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv")

P25_file <- file.path(prior_dir,
                      "HER2_UP_overlap_25_FDR0.05_LFC1_ALL3_with_stats_Wolf_Robinson_Brueffer.csv")

# --------- Load Expression ---------
Xdf <- fread(expr_file)
genes <- Xdf[[1]]
Xdf[[1]] <- NULL

expr_mat <- as.matrix(Xdf)
rownames(expr_mat) <- genes
storage.mode(expr_mat) <- "numeric"

# --------- Z-score within cohort ---------
Ez <- t(scale(t(expr_mat)))

# --------- Load Gene Sets ---------
P76_genes <- fread(P76_file)$gene_symbol
P25_genes <- fread(P25_file)$gene_symbol

P76_genes <- intersect(P76_genes, rownames(Ez))
P25_genes <- intersect(P25_genes, rownames(Ez))

cat("Matched P76 genes:", length(P76_genes), "\n")
cat("Matched P25 genes:", length(P25_genes), "\n")

# --------- Compute Scores ---------
P76_meanZ <- colMeans(Ez[P76_genes, , drop = FALSE], na.rm = TRUE)
P25_meanZ <- colMeans(Ez[P25_genes, , drop = FALSE], na.rm = TRUE)

score_df <- data.frame(
  SampleID = colnames(Ez),
  P76_meanZ = P76_meanZ,
  P25_meanZ = P25_meanZ
)

# --------- Merge with HER2 Status ---------
meta_used <- fread(meta_file)
score_df <- merge(score_df, meta_used[, .(SampleID, Group)], by = "SampleID")

# --------- Save Scores ---------
out_scores <- file.path(res_dir, "Robinson_P76P25_scores_meanZ.tsv")
fwrite(score_df, out_scores, sep = "\t")
cat("Saved:", out_scores, "\n")

# --------- Boxplot Validation ---------
p <- ggplot(score_df, aes(x = Group, y = P76_meanZ, fill = Group)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "Robinson: P76 meanZ by HER2 status",
       y = "P76 meanZ score")

ggsave(file.path(res_dir, "Robinson_P76_meanZ_boxplot.png"),
       plot = p, width = 6, height = 5, dpi = 220)

# --------- Wilcoxon Test ---------
w <- wilcox.test(P76_meanZ ~ Group, data = score_df)
cat("Wilcoxon p-value (P76):", signif(w$p.value, 4), "\n")

# --------- ERBB2 Correlation ---------
if ("ERBB2" %in% rownames(expr_mat)) {
  erbb2 <- expr_mat["ERBB2", ]
  rho <- cor(P76_meanZ, erbb2, method = "spearman",
             use = "pairwise.complete.obs")
  cat("Spearman rho with ERBB2:", signif(rho, 4), "\n")
}

cat("Scoring complete.\n")
