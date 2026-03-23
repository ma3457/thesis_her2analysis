# ================= run_brueffer_joint_scores.R =================
# Merge Brueffer P76/P25 + HARPS per-sample scores and plot P76 vs HARPS.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

base_dir <- "~/Desktop/Thesis/AIM 2 HER2 DATA/brueffer et al"
out_dir  <- file.path(base_dir, "aim2_outputs", "GSE81538")

p76p25_fp <- file.path(out_dir, "Brueffer_P76P25_scores_meanZ.tsv")
harps_fp  <- file.path(out_dir, "Brueffer_HARPS_scores_meanZ.tsv")

stopifnot(file.exists(p76p25_fp), file.exists(harps_fp))

p76p25 <- fread(p76p25_fp)
harps  <- fread(harps_fp)

p76p25$SampleID <- as.character(p76p25$SampleID)
harps$SampleID  <- as.character(harps$SampleID)

joint <- left_join(p76p25, harps[, .(SampleID, HARPS_meanZ)], by = "SampleID") %>% as.data.table()

out_fp <- file.path(out_dir, "Brueffer_joint_scores_P76P25_HARPS.tsv")
fwrite(joint, out_fp, sep = "\t")
message("Saved: ", out_fp)

p <- ggplot(joint, aes(x = P76_meanZ, y = HARPS_meanZ, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_classic(base_size = 18) +
  labs(
    title = "Brueffer: P76 vs HARPS (two-axis HER2 signaling)",
    x = "P76 meanZ (amplification axis)",
    y = "HARPS meanZ (HER2-like axis)"
  )

ggsave(file.path(out_dir, "Brueffer_P76_vs_HARPS_scatter.png"),
       p, width = 10, height = 7, dpi = 220)

message("DONE: Brueffer joint scoring complete.")
# ================= end =================
