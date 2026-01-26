# ================= Debets HER2/HARP validation barplot =================
# This makes a simple PDF barplot summarizing:
#  - total prioritized HER2/HARP features
#  - detected in Debets proteome
#  - nominally significant (P < 0.05)
#  - FDR-significant (FDR < 0.05)

base_dir <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/Debets et al/aim2_outputs/Debets_HER2_HARP_validation"
out_pdf  <- file.path(base_dir, "Debets_HER2_HARP_validation_barplot.pdf")

# Manually set the counts based on the validation summary
counts <- c(
  total_prioritized   = 12,  # from Wolf core list
  detected_in_Debets  = 8,   # present in Debets proteome
  nominal_P_lt_0.05   = 1,   # nominal hit(s)
  FDR_sig_lt_0.05     = 0    # none pass FDR
)

labels <- c("Prioritized\nfeatures",
            "Detected in\nDebets proteome",
            "Nominally\nsignificant (P<0.05)",
            "FDR-significant\n(FDR<0.05)")

pdf(out_pdf, width = 7, height = 5)
par(mar = c(7, 5, 4, 2))  # extra bottom margin for wrapped labels

bp <- barplot(
  counts,
  names.arg = labels,
  ylim = c(0, max(counts) + 2),
  ylab = "Number of features",
  main = "Debets – HER2/HARP Validation Summary",
  cex.names = 0.9,
  cex.axis = 0.9
)

# Add value labels on top of bars
text(
  x = bp,
  y = counts,
  labels = counts,
  pos = 3,
  cex = 0.9
)

dev.off()

message("Wrote barplot → ", out_pdf)
# ================= END =================
