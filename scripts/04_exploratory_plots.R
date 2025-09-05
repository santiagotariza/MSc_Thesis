###############################################################################
# Script: 04_exploratory_plots.R
# Purpose: Generate exploratory plots from variance-stabilized data (VST):
#   - PCA plots
#   - Heatmaps of top variable genes
#
# Dependencies:
#   - DESeq2
#   - ggplot2
#   - pheatmap
#   - matrixStats
#   - ggrepel
#   - dplyr
#   - SummarizedExperiment
#   - rstudioapi
#
# Inputs:
#   - results/01_filtered_data.RData (sample_info, counts_mRNA, counts_lncRNA)
#   - results/03_dds_vsd_objects.RData (mRNA_obj, lncRNA_obj)
#
# Outputs:
#   - results/04_mRNA_pca_vsd.png
#   - results/04_mRNA_heatmap_top_var_genes_vsd.png
#   - results/04_lncRNA_pca_vsd.png
#   - results/04_lncRNA_heatmap_top_var_genes_vsd.png
###############################################################################

# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(matrixStats)
library(ggrepel)
library(dplyr)
library(SummarizedExperiment)
library(rstudioapi)

# Define directories
script_dir  <- dirname(rstudioapi::getSourceEditorContext()$path)
results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Load data from previous steps
# -------------------------------------------------------------------------
load(file.path(results_dir, "01_filtered_data.RData"))
load(file.path(results_dir, "03_dds_vsd_objects.RData"))

# -------------------------------------------------------------------------
# Function: PCA and Heatmap plots
# -------------------------------------------------------------------------
plot_pca_and_heatmap <- function(vsd_obj, sample_info, prefix) {
  if (is.null(vsd_obj) || nrow(assay(vsd_obj)) < 2) {
    warning("Skipping PCA/Heatmap for ", prefix, 
            ": VSD object is empty or has too few genes.")
    return(NULL)
  }
  
  message("Generating PCA and Heatmap for ", prefix, " (uncorrected)...")
  
  # Ensure sample_info and vsd_obj are aligned
  sample_info_matched <- sample_info[colnames(assay(vsd_obj)), , drop = FALSE]
  
  # --- PCA Plot ---
  pca <- prcomp(t(assay(vsd_obj)))
  pca_data <- as.data.frame(pca$x)
  pca_data$condition   <- sample_info_matched$condition
  pca_data$sample_name <- rownames(sample_info_matched)
  
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
  
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = sample_name), size = 3, box.padding = 0.5) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(paste0("PCA of ", prefix, " (uncorrected)")) +
    theme_minimal() +
    theme(legend.position = "right")
  
  ggsave(file.path(results_dir, paste0("04_", prefix, "_pca_vsd.png")), 
         plot = pca_plot, width = 8, height = 7, dpi = 300)
  
  message("Saved PCA plot for ", prefix)
  
  # --- Heatmap Plot ---
  num_genes_for_heatmap <- min(50, nrow(assay(vsd_obj)))
  if (num_genes_for_heatmap < 2) {
    warning("Less than 2 variable genes available for heatmap of ", prefix, 
            ". Skipping heatmap.")
    return(NULL)
  }
  
  top_var_genes <- head(order(rowVars(assay(vsd_obj)), decreasing = TRUE), 
                        num_genes_for_heatmap)
  mat <- assay(vsd_obj)[top_var_genes, ]
  
  annotation_col <- as.data.frame(colData(vsd_obj)[, "condition", drop = FALSE])
  rownames(annotation_col) <- colnames(mat)
  
  png(file.path(results_dir, paste0("04_", prefix, "_heatmap_top_var_genes_vsd.png")), 
      width = 1200, height = 1200, res = 150)
  pheatmap(
    mat,
    cluster_rows   = TRUE,
    show_rownames  = FALSE,
    cluster_cols   = TRUE,
    annotation_col = annotation_col,
    main           = paste0("Heatmap of top ", num_genes_for_heatmap, 
                            " variable genes (", prefix, ", VST)"),
    fontsize       = 8
  )
  dev.off()
  
  message("Saved heatmap for ", prefix)
}

# -------------------------------------------------------------------------
# Run for mRNA
# -------------------------------------------------------------------------
if (!is.null(mRNA_obj) && "vsd" %in% names(mRNA_obj)) {
  plot_pca_and_heatmap(mRNA_obj$vsd, sample_info, "mRNA")
} else {
  warning("mRNA_obj$vsd is not available. Skipping mRNA analysis.")
}

# -------------------------------------------------------------------------
# Run for lncRNA
# -------------------------------------------------------------------------
if (!is.null(lncRNA_obj) && "vsd" %in% names(lncRNA_obj)) {
  plot_pca_and_heatmap(lncRNA_obj$vsd, sample_info, "lncRNA")
} else {
  warning("lncRNA_obj$vsd is not available. Skipping lncRNA analysis.")
}

# -------------------------------------------------------------------------
# End
# -------------------------------------------------------------------------
message("Step 04 completed. Exploratory plots generated. Check the 'results' folder for PNG outputs.")