###############################################################################
# Script: 02_exploratory_qc_analysis.R
# Purpose: Perform exploratory quality control (QC) analysis, including:
#   - Parse and visualize Fastp metrics (Q20/Q30 rates).
#   - Plot raw counts distributions (boxplot + density).
#   - Generate correlation heatmaps of samples.
#
# Dependencies:
#   - tidyverse
#   - jsonlite
#   - pheatmap
#   - RColorBrewer
#   - cowplot
#   - rstudioapi
#
# Inputs:
#   - results/01_filtered_data.RData (sample_info, counts_mRNA, counts_lncRNA)
#   - fastp_reports/*.json (Fastp summary reports)
#
# Outputs:
#   - results/02_fastp_quality_metrics.png
#   - results/02_mRNA_raw_counts_distribution.png
#   - results/02_lncRNA_raw_counts_distribution.png
#   - results/02_mRNA_correlation_heatmap.png
#   - results/02_lncRNA_correlation_heatmap.png
###############################################################################

# Load required libraries
library(tidyverse)
library(jsonlite)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(rstudioapi)

# Define directories
script_dir  <- dirname(rstudioapi::getSourceEditorContext()$path)
results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Directory containing Fastp JSON reports
fastp_json_dir <- file.path(script_dir, "../fastp_reports")

# -------------------------------------------------------------------------
# Diagnostics
# -------------------------------------------------------------------------
message("Searching for Fastp JSON reports in: ", fastp_json_dir)
json_files <- list.files(fastp_json_dir, pattern = "\\.json$", full.names = TRUE)
message("Number of Fastp JSON files found: ", length(json_files))

if (length(json_files) == 0) {
  warning("No Fastp JSON files found. Skipping Fastp metrics visualization.")
}

# -------------------------------------------------------------------------
# 1. Load filtered data (from previous script)
# -------------------------------------------------------------------------
load(file.path(results_dir, "01_filtered_data.RData"))

# -------------------------------------------------------------------------
# 2. Fastp metrics analysis
# -------------------------------------------------------------------------
message("Analyzing Fastp reports...")
json_files <- list.files(fastp_json_dir, pattern = "\\.json$", full.names = TRUE)

if (length(json_files) > 0) {
  fastp_metrics <- map_dfr(json_files, function(file) {
    data <- fromJSON(file)
    tibble(
      sample          = sub(".fastp.json", "", basename(file)),
      q20_rate_before = data$summary$before_filtering$q20_rate,
      q30_rate_before = data$summary$before_filtering$q30_rate,
      q20_rate_after  = data$summary$after_filtering$q20_rate,
      q30_rate_after  = data$summary$after_filtering$q30_rate
    )
  })
  
  # Align metrics with sample_info
  fastp_metrics_matched <- left_join(sample_info, fastp_metrics, by = "sample")
  
  # Barplot of Q30 after filtering
  p <- ggplot(fastp_metrics_matched, aes(x = reorder(sample, q30_rate_after), 
                                         y = q30_rate_after, fill = condition)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.2f%%", q30_rate_after * 100)), 
              vjust = -0.5, size = 3) +
    labs(
      title = "Q30 base rate after filtering",
      x = "Sample",
      y = "Q30 rate"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave(file.path(results_dir, "02_fastp_quality_metrics.png"), 
         plot = p, width = 10, height = 7)
  
  message("Saved Fastp metrics plot to results/02_fastp_quality_metrics.png")
}

# -------------------------------------------------------------------------
# 3. Raw counts distribution plots
# -------------------------------------------------------------------------
plot_raw_counts_distribution <- function(counts_matrix, results_dir, type_name) {
  counts_df <- as.data.frame(t(counts_matrix))
  
  counts_long <- counts_df %>%
    mutate(sample = rownames(counts_df)) %>%
    pivot_longer(cols = -sample, names_to = "gene", values_to = "count") %>%
    mutate(count_log2 = log2(count + 1))
  
  # Boxplot
  p_boxplot <- ggplot(counts_long, aes(x = sample, y = count_log2, fill = sample)) +
    geom_boxplot() +
    labs(
      title = paste("Gene counts distribution (", type_name, ")", sep = ""),
      x = "Sample",
      y = "Count (log2)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          legend.position = "none")
  
  # Density plot
  p_density <- ggplot(counts_long, aes(x = count_log2, color = sample)) +
    geom_density(alpha = 0.5, linewidth = 1) +
    labs(
      title = paste("Gene counts density distribution (", type_name, ")", sep = ""),
      x = "Count (log2)",
      y = "Density"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Combine plots
  combined_plot <- cowplot::plot_grid(p_boxplot, p_density, nrow = 1)
  
  ggsave(file.path(results_dir, paste0("02_", type_name, "_raw_counts_distribution.png")), 
         plot = combined_plot, width = 12, height = 6)
  
  message("Saved raw counts distribution plots for ", type_name)
}

# Generate for mRNA and lncRNA
plot_raw_counts_distribution(counts_mRNA, results_dir, "mRNA")
plot_raw_counts_distribution(counts_lncRNA, results_dir, "lncRNA")

# -------------------------------------------------------------------------
# 4. Correlation heatmap
# -------------------------------------------------------------------------
plot_correlation_heatmap <- function(counts_matrix, sample_info, type_name, results_dir) {
  counts_matrix <- counts_matrix[, colnames(counts_matrix) %in% rownames(sample_info)]
  sample_info_matched <- sample_info[colnames(counts_matrix), ]
  
  mat_log2 <- log2(counts_matrix + 1)
  cor_mat  <- cor(mat_log2)
  
  annotation_col <- dplyr::select(sample_info_matched, condition)
  rownames(annotation_col) <- sample_info_matched$sample
  
  png(file.path(results_dir, paste0("02_", type_name, "_correlation_heatmap.png")),
      width = 1200, height = 1200, res = 150)
  
  pheatmap(
    cor_mat,
    main            = paste("Sample correlation heatmap (", type_name, ", log2(counts+1))", sep = ""),
    annotation_col  = annotation_col,
    fontsize_row    = 8,
    fontsize_col    = 8,
    color           = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
    display_numbers = TRUE,
    fontsize_number = 6,
    number_color    = "black"
  )
  
  dev.off()
  
  message("Saved correlation heatmap for ", type_name)
}

# Generate for mRNA and lncRNA
plot_correlation_heatmap(counts_mRNA, sample_info, "mRNA", results_dir)
plot_correlation_heatmap(counts_lncRNA, sample_info, "lncRNA", results_dir)

# -------------------------------------------------------------------------
# End
# -------------------------------------------------------------------------
message("Step 02 completed. QC analysis finished. Check the 'results' folder for PNG outputs.")