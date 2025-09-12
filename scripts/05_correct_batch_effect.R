##############################################
# Script: 05_correct_batch_effect.R
# Purpose: Apply batch effect correction (ComBat)
# Input:   01_filtered_data.RData, 03_dds_vsd_objects.RData
# Output:  PCA plots, heatmaps, 05_combat_corrected_vsd_objects.RData
##############################################

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("sva", update = TRUE, force =TRUE)


# --- Libraries ---
library(tidyverse)
library(DESeq2)
library(sva)                # For ComBat
library(ggplot2)
library(pheatmap)
library(matrixStats)
library(ggrepel)
library(SummarizedExperiment) # For assay()

# --- Define directories ---
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load data ---
load(file.path(results_dir, "01_filtered_data.RData"))
load(file.path(results_dir, "03_dds_vsd_objects.RData"))

# --- Define batch variable from sample names ---
message("\n--- Defining batch variable from sample IDs ---")
sample_info$batch <- factor(sub(".*([0-9]+)$", "\\1", rownames(sample_info)))

# Ensure sample_info is aligned with counts
if (!all(rownames(sample_info) %in% colnames(counts_mRNA))) {
  warning("Sample names in 'sample_info' do not match 'counts_mRNA'. Check '01_load_filter_data.R'. Aligning automatically.")
  sample_info <- sample_info[colnames(counts_mRNA), ]
}

# --- Function: Apply ComBat correction and generate plots ---
process_and_plot_corrected_data <- function(vsd_obj, sample_info_df, prefix) {
  if (is.null(vsd_obj) || nrow(assay(vsd_obj)) < 2) {
    message(paste("Skipping", prefix, ": empty or insufficient genes in VSD object."))
    return(NULL)
  }
  
  message(paste0("\n--- Processing ", prefix, " for batch correction ---"))
  
  # Align sample info with VSD
  sample_info_df_matched <- sample_info_df[colnames(assay(vsd_obj)), , drop = FALSE]
  mat_vsd <- assay(vsd_obj)
  
  # --- Remove genes with near-zero variance ---
  row_vars <- rowVars(mat_vsd)
  genes_to_keep <- which(row_vars > 1e-6)
  if (length(genes_to_keep) == 0) {
    message(paste("No genes with sufficient variance for", prefix, ". Skipping PCA/heatmap."))
    return(NULL)
  }
  mat_vsd_filtered <- mat_vsd[genes_to_keep, , drop = FALSE]
  
  # --- Apply ComBat ---
  message(paste("Applying ComBat to", prefix, "VSD data..."))
  mod_combat <- tryCatch({
    model.matrix(~ condition, data = sample_info_df_matched)
  }, error = function(e) {
    message("  Error in model.matrix(~condition). Falling back to ~1.")
    model.matrix(~1, data = sample_info_df_matched)
  })
  
  # Check model matrix rank
  if (qr(mod_combat)$rank < ncol(mod_combat)) {
    message("  Warning: Model matrix not full rank. Falling back to ~1.")
    mod_combat <- model.matrix(~1, data = sample_info_df_matched)
  }
  
  mat_corrected <- ComBat(dat = mat_vsd_filtered,
                          batch = sample_info_df_matched$batch,
                          mod = mod_combat)
  
  # --- Remove problematic rows after ComBat ---
  if (any(is.na(mat_corrected)) || any(is.infinite(mat_corrected))) {
    mat_corrected <- mat_corrected[!apply(is.na(mat_corrected) | is.infinite(mat_corrected), 1, any), , drop = FALSE]
  }
  
  if (nrow(mat_corrected) < 2 || ncol(mat_corrected) < 2) {
    message(paste("Corrected matrix for", prefix, "has <2 genes or samples. Skipping PCA."))
    return(NULL)
  }
  
  # --- PCA plot ---
  pca <- prcomp(t(mat_corrected))
  pca_data <- as.data.frame(pca$x)
  samples_remaining <- colnames(mat_corrected)
  sample_info_df_pca <- sample_info_df_matched[samples_remaining, , drop = FALSE]
  
  pca_data$condition <- sample_info_df_pca$condition
  pca_data$batch <- sample_info_df_pca$batch
  pca_data$sample_name <- rownames(sample_info_df_pca)
  
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
  
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, shape = batch)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = sample_name), size = 3, box.padding = 0.5) +
    xlab(paste0("PC1 (", percentVar[1], "%)")) +
    ylab(paste0("PC2 (", percentVar[2], "%)")) +
    ggtitle(paste0(prefix, " PCA (ComBat corrected)")) +
    theme_minimal() +
    theme(legend.position = "right")
  
  ggsave(file.path(results_dir, paste0("05_", prefix, "_pca_combat_corrected.png")),
         plot = pca_plot, width = 8, height = 7, dpi = 300)
  
  # --- Heatmap plot ---
  num_genes <- min(50, nrow(mat_corrected))
  if (num_genes >= 2) {
    top_var_genes <- head(order(rowVars(mat_corrected), decreasing = TRUE), num_genes)
    mat_for_heatmap <- mat_corrected[top_var_genes, , drop = FALSE]
    annotation_col <- as.data.frame(sample_info_df_pca[, c("condition", "batch"), drop = FALSE])
    
    png(file.path(results_dir, paste0("05_", prefix, "_heatmap_combat_corrected.png")),
        width = 1200, height = 1200, res = 150)
    pheatmap(mat_for_heatmap,
             cluster_rows = TRUE,
             show_rownames = FALSE,
             cluster_cols = TRUE,
             annotation_col = annotation_col,
             main = paste0("Top ", num_genes, " variable genes (", prefix, ", ComBat corrected)"),
             fontsize = 8)
    dev.off()
  }
  
  return(mat_corrected)
}

# --- Run for mRNA ---
if (!is.null(mRNA_obj) && "vsd" %in% names(mRNA_obj)) {
  corrected_matrix_mRNA <- process_and_plot_corrected_data(mRNA_obj$vsd, sample_info, "mRNA")
} else {
  message("Skipping mRNA: 'mRNA_obj$vsd' not available.")
  corrected_matrix_mRNA <- NULL
}

# --- Run for lncRNA ---
if (!is.null(lncRNA_obj) && "vsd" %in% names(lncRNA_obj)) {
  corrected_matrix_lncRNA <- process_and_plot_corrected_data(lncRNA_obj$vsd, sample_info, "lncRNA")
} else {
  message("Skipping lncRNA: 'lncRNA_obj$vsd' not available.")
  corrected_matrix_lncRNA <- NULL
}

# --- Save corrected VSD objects ---
if (!is.null(corrected_matrix_mRNA)) {
  vsd_corrected_mRNA <- mRNA_obj$vsd
  assay(vsd_corrected_mRNA) <- corrected_matrix_mRNA
} else {
  vsd_corrected_mRNA <- NULL
}

if (!is.null(corrected_matrix_lncRNA)) {
  vsd_corrected_lncRNA <- lncRNA_obj$vsd
  assay(vsd_corrected_lncRNA) <- corrected_matrix_lncRNA
} else {
  vsd_corrected_lncRNA <- NULL
}

save(vsd_corrected_mRNA, vsd_corrected_lncRNA,
     file = file.path(results_dir, "05_combat_corrected_vsd_objects.RData"))

message("\n--- Batch correction completed. Results saved to '05_combat_corrected_vsd_objects.RData' ---")