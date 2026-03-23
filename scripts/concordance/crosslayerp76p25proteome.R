# ================== crosslayerp76p25proteome.R ==================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------

BASE_DIR <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA"

OUT_DIR  <- file.path(BASE_DIR, "crosslayerP76P25validation")
TAB_DIR  <- file.path(OUT_DIR, "tables")
PLOT_DIR <- file.path(OUT_DIR, "plots")
AUD_DIR  <- file.path(OUT_DIR, "audits")

dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(AUD_DIR,  showWarnings = FALSE, recursive = TRUE)

# Inputs
P76_FILE <- file.path(BASE_DIR, "aim2_prioritized",
                      "HER2_UP_overlap_76_with_stats_Wolf_Robinson_Brueffer.csv")

P25_FILE <- file.path(BASE_DIR, "aim2_prioritized",
                      "HER2_UP_overlap_25_FDR0.05_LFC1_ALL3_with_stats_Wolf_Robinson_Brueffer.csv")

RNA_DE_FILE <- file.path(BASE_DIR, "wolf et al", "aim2_outputs",
                         "GSE194040_ISPY2_mRNA_BPsubtype",
                         "ALL_HER2pos_vs_HER2neg_DE.tsv")

RAJKUMAR_DE_FILE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/raj kumar et al/aim2_outputs/RajKumar_PROT_GLOBAL/DE_full.tsv"

KRUG_DE_FILE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/krug et al/aim2_outputs/Krug_CPTAC_BRCA_PROT_A_ClinicalOnly/DE_full.tsv"

MERTINS_DE_FILE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/mertins et al/aim2_outputs/Mertins_CPTAC_BRCA_PROT_GLOBAL_ClinicalOnly/DE_full.tsv"

# ------------------------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------------------------

AMP_GENES <- c("ERBB2","GRB7","STARD3","MIEN1","PGAP3","ORMDL3","PSMD3","IKZF3")

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------

clean_gene <- function(x){
  x <- as.character(x)
  x <- str_trim(x)
  x <- gsub("\\.\\d+$","",x)
  toupper(x)
}

pick_first_present <- function(nms, cands){
  hit <- intersect(cands,nms)
  if(length(hit)) return(hit[1])
  NA_character_
}

read_gene_list <- function(path){
  dt <- fread(path)
  g <- pick_first_present(names(dt), c("gene_symbol","feature_id","Gene","gene","SYMBOL","symbol"))
  if (is.na(g)) stop("No gene column found in: ", path)
  unique(clean_gene(dt[[g]]))
}

read_de_standard <- function(path, label){
  dt <- fread(path)
  
  gene_col <- pick_first_present(
    names(dt),
    c("gene_symbol","feature_id","Gene","gene","SYMBOL","symbol","Protein","protein")
  )
  
  logfc_col <- pick_first_present(
    names(dt),
    c("logFC","log2FC","log2FoldChange","log2FC_HER2pos_vs_HER2neg","logFC_HER2pos_vs_HER2neg")
  )
  
  p_col <- pick_first_present(
    names(dt),
    c("P.Value","PValue","p.value","pval","PValue_HER2pos_vs_HER2neg","pvalue_HER2pos_vs_HER2neg","P.Value_HER2pos_vs_HER2neg")
  )
  
  fdr_col <- pick_first_present(
    names(dt),
    c("adj.P.Val","adj.P.Val.","padj","FDR","qval","adj_pval","FDR_HER2pos_vs_HER2neg","adj.P.Val_HER2pos_vs_HER2neg","padj_HER2pos_vs_HER2neg")
  )
  
  if (is.na(gene_col)) stop("No gene column found in: ", path)
  if (is.na(logfc_col)) stop("No logFC column found in: ", path)
  if (is.na(p_col)) stop("No P.Value column found in: ", path)
  if (is.na(fdr_col)) stop("No adj.P.Val/FDR column found in: ", path)
  
  out <- data.table(
    gene_symbol = clean_gene(dt[[gene_col]]),
    logFC = suppressWarnings(as.numeric(dt[[logfc_col]])),
    P.Value = suppressWarnings(as.numeric(dt[[p_col]])),
    adj.P.Val = suppressWarnings(as.numeric(dt[[fdr_col]]))
  )
  
  out[, abs_logFC := abs(logFC)]
  setorder(out, gene_symbol, adj.P.Val, P.Value, -abs_logFC)
  out <- out[!duplicated(gene_symbol)]
  out[, abs_logFC := NULL]
  
  setnames(
    out,
    old = c("logFC","P.Value","adj.P.Val"),
    new = paste0(label, c("_logFC","_P.Value","_adj.P.Val"))
  )
  
  out
}

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------

P76 <- read_gene_list(P76_FILE)
P25 <- read_gene_list(P25_FILE)

rna  <- read_de_standard(RNA_DE_FILE, "RNA")
raj  <- read_de_standard(RAJKUMAR_DE_FILE, "Raj")
krug <- read_de_standard(KRUG_DE_FILE, "Krug")
mert <- read_de_standard(MERTINS_DE_FILE, "Mert")

# ------------------------------------------------------------------------------
# CORE FUNCTION
# ------------------------------------------------------------------------------

run_analysis <- function(genes, name, prot, label){
  
  df <- data.table(gene_symbol = genes)
  df <- merge(df, rna,  by = "gene_symbol", all.x = TRUE)
  df <- merge(df, prot, by = "gene_symbol", all.x = TRUE)
  
  prot_col <- paste0(label, "_logFC")
  
  df[, Concordance := ifelse(
    is.na(get(prot_col)), "Not_detected",
    ifelse(RNA_logFC > 0 & get(prot_col) > 0, "Concordant",
           ifelse(RNA_logFC < 0 & get(prot_col) < 0, "Concordant", "Discordant"))
  )]
  
  df[, Amplicon := ifelse(gene_symbol %in% AMP_GENES, "Amplicon", "Non-amplicon")]
  
  # Save merged table
  out_csv <- file.path(TAB_DIR, paste0(name, "_", label, ".csv"))
  fwrite(df, out_csv)
  
  # Console summary
  cat("\n", name, " - ", label, "\n", sep = "")
  print(table(df$Concordance, useNA = "ifany"))
  print(table(df$Amplicon, df$Concordance, useNA = "ifany"))
  
  # Save audit text
  audit_file <- file.path(AUD_DIR, paste0(name, "_", label, "_audit.txt"))
  audit_lines <- capture.output({
    cat(name, "-", label, "\n")
    print(table(df$Concordance, useNA = "ifany"))
    print(table(df$Amplicon, df$Concordance, useNA = "ifany"))
  })
  writeLines(audit_lines, audit_file)
  
  # Plot dataframe
  plot_df <- df[!is.na(df[[prot_col]]), ]
  
  if (nrow(plot_df) > 0) {
    p <- ggplot(
      plot_df,
      aes(
        x = RNA_logFC,
        y = .data[[prot_col]],
        color = Concordance,
        shape = Amplicon
      )
    ) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(
        title = paste0(name, " RNA vs Protein (", label, ")"),
        x = "RNA log2FC",
        y = "Protein log2FC"
      ) +
      theme_classic()
    
    ggsave(
      file.path(PLOT_DIR, paste0(name, "_", label, "_scatter.png")),
      p,
      width = 7,
      height = 5,
      dpi = 300
    )
    
    # Optional concordance bar chart
    summary_df <- df %>%
      group_by(Concordance) %>%
      summarise(n = n(), .groups = "drop")
    
    p_bar <- ggplot(summary_df, aes(x = Concordance, y = n, fill = Concordance)) +
      geom_col() +
      labs(
        title = paste0(name, " Concordance Summary (", label, ")"),
        x = NULL,
        y = "Gene count"
      ) +
      theme_classic()
    
    ggsave(
      file.path(PLOT_DIR, paste0(name, "_", label, "_bar.png")),
      p_bar,
      width = 6,
      height = 4,
      dpi = 300
    )
  } else {
    cat("No proteins detected for plotting in ", name, " / ", label, "\n", sep = "")
  }
}

# ------------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------------

run_analysis(P76, "P76", raj,  "Raj")
run_analysis(P76, "P76", krug, "Krug")
run_analysis(P76, "P76", mert, "Mert")

run_analysis(P25, "P25", raj,  "Raj")
run_analysis(P25, "P25", krug, "Krug")
run_analysis(P25, "P25", mert, "Mert")

cat("\nDONE\n")
# =====================================================================
# ---- PLOTTING FROM SAVED TABLES (APPENDED SECTION) ----

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

CROSSLAYER_BASE <- "/Users/maya.anand/Desktop/Thesis/AIM 2 HER2 DATA/crosslayerP76P25validation"
TABLE_DIR <- file.path(CROSSLAYER_BASE, "tables")
PLOT_DIR  <- file.path(CROSSLAYER_BASE, "plots")

make_crosslayer_plots <- function(file_path, dataset_label, program_label) {
  
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # identify protein column
  prot_col <- grep("_logFC$", names(df), value = TRUE)
  prot_col <- prot_col[prot_col != "RNA_logFC"][1]
  
  df <- df %>%
    mutate(
      Concordance = factor(Concordance, levels = c("Concordant", "Discordant", "Not_detected")),
      Amplicon = factor(Amplicon, levels = c("Amplicon", "Non-amplicon"))
    )
  
  # ---- BAR ----
  p_bar <- df %>%
    count(Concordance) %>%
    ggplot(aes(x = Concordance, y = n, fill = Concordance)) +
    geom_col() +
    theme_classic() +
    labs(
      title = paste0(program_label, " Concordance Summary (", dataset_label, ")"),
      y = "Gene count"
    )
  
  ggsave(file.path(PLOT_DIR, paste0(program_label, "_", dataset_label, "_bar.png")),
         p_bar, width = 7, height = 5)
  
  # ---- SCATTER ----
  scatter_df <- df %>%
    filter(!is.na(RNA_logFC), !is.na(.data[[prot_col]]))
  
  if (nrow(scatter_df) > 0) {
    p_scatter <- ggplot(
      scatter_df,
      aes(
        x = RNA_logFC,
        y = .data[[prot_col]],
        color = Concordance,
        shape = Amplicon
      )
    ) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_classic() +
      labs(
        title = paste0(program_label, " RNA vs Protein (", dataset_label, ")"),
        x = "RNA log2FC",
        y = "Protein log2FC"
      )
    
    ggsave(file.path(PLOT_DIR, paste0(program_label, "_", dataset_label, "_scatter.png")),
           p_scatter, width = 7, height = 5)
  }
  
  cat("Plotted:", program_label, dataset_label, "\n")
}

# ---- RUN PLOTS ----
make_crosslayer_plots(file.path(TABLE_DIR, "P25_Raj.csv"),  "Raj",  "P25")
make_crosslayer_plots(file.path(TABLE_DIR, "P25_Krug.csv"), "Krug", "P25")
make_crosslayer_plots(file.path(TABLE_DIR, "P25_Mert.csv"), "Mert", "P25")

make_crosslayer_plots(file.path(TABLE_DIR, "P76_Raj.csv"),  "Raj",  "P76")
make_crosslayer_plots(file.path(TABLE_DIR, "P76_Krug.csv"), "Krug", "P76")
make_crosslayer_plots(file.path(TABLE_DIR, "P76_Mert.csv"), "Mert", "P76")

cat("\nALL PLOTS COMPLETE\n")