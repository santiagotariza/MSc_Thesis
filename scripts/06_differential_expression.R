##############################################
# Script: 06_differential_expression.R
# Purpose: Perform differential expression analysis with DESeq2
# Input:   03_dds_vsd_objects.RData, 05_combat_corrected_vsd_objects.RData
# Output:  Differential expression results (CSV), volcano plots (PNG)
##############################################

# --- Libraries ---
library(DESeq2)
library(apeglm)
library(dplyr)
library(ggplot2)
library(ggrepel)

# --- Parameters ---
log2FC_threshold <- 1
padj_threshold <- 0.01
use_lfc_shrink <- TRUE
volcano_xlim <- c(-7, 7)
volcano_ylim_max <- 25
top_genes_to_label <- 10

# --- Define directories ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}
results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load data ---
load(file.path(results_dir, "03_dds_vsd_objects.RData"))
load(file.path(results_dir, "05_combat_corrected_vsd_objects.RData"))

# --- Define comparisons of interest ---
comparisons <- list(
  list(type = "name", numerator = "G3_HMGB1_KO", denominator = "G1_NTC"),
  list(type = "name", numerator = "G4_HMGB2_KO", denominator = "G1_NTC"),
  list(type = "name", numerator = "G2_WTC", denominator = "G1_NTC"),
  list(type = "contrast", numerator = "G3_HMGB1_KO", denominator = "G4_HMGB2_KO"),
  list(type = "contrast", numerator = "G3_HMGB1_KO", denominator = "G2_WTC"),
  list(type = "contrast", numerator = "G4_HMGB2_KO", denominator = "G2_WTC")
)

# --- Function: Run DESeq2 analysis and generate volcano plots ---
run_deseq_analysis <- function(dds_obj, vsd_corrected_obj, biotype) {
  if (is.null(dds_obj) || !("dds" %in% names(dds_obj))) {
    message(paste("Skipping", biotype, ": DESeq2 object not available."))
    return(NULL)
  }
  
  dds <- dds_obj$dds
  message(paste("\n--- Running differential expression analysis for", biotype, "---"))
  
  for (comp in comparisons) {
    numerator <- comp$numerator
    denominator <- comp$denominator
    comp_name <- paste0(numerator, "_vs_", denominator)
    
    message(paste("  - Comparison:", comp_name))
    
    # Run DESeq2 results
    if (comp$type == "name") {
      res <- results(dds, name = paste0("condition_", numerator, "_vs_", denominator), alpha = padj_threshold)
    } else {
      res <- results(dds, contrast = c("condition", numerator, denominator), alpha = padj_threshold)
    }
    
    # Apply LFC shrinkage
    if (use_lfc_shrink) {
      if (comp$type == "name") {
        message("  - Applying lfcShrink with 'apeglm'...")
        res <- lfcShrink(dds, coef = paste0("condition_", numerator, "_vs_", denominator), type = "apeglm")
      } else {
        message("  - Applying lfcShrink with 'normal' (fallback)...")
        res <- lfcShrink(dds, contrast = c("condition", numerator, denominator), type = "normal")
        warning("'apeglm' not supported for contrast. Using 'normal'.")
      }
    }
    
    # Filter significant results
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      filter(!is.na(padj) & abs(log2FoldChange) >= log2FC_threshold & padj < padj_threshold)
    
    output_filename <- file.path(results_dir, paste0("06_", biotype, "_", comp_name, "_DE_results.csv"))
    write.csv(res_df, file = output_filename, row.names = FALSE)
    
    message(paste("  - Found", nrow(res_df), "significantly DE genes/lncRNAs"))
    message(paste("  - Results saved to:", output_filename))
    
    # Volcano plot (using corrected VSD object only for visualization consistency)
    res_plot <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      mutate(sig = case_when(
        log2FoldChange > log2FC_threshold & padj < padj_threshold ~ "Up-regulated",
        log2FoldChange < -log2FC_threshold & padj < padj_threshold ~ "Down-regulated",
        TRUE ~ "Not significant"
      ))
    
    res_labeled <- res_plot %>%
      filter(sig != "Not significant") %>%
      arrange(padj) %>%
      head(top_genes_to_label)
    
    volcano_plot <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not significant" = "grey")) +
      geom_text_repel(
        data = res_labeled,
        aes(label = gene_id),
        size = 3,
        box.padding = unit(0.5, "lines"),
        point.padding = unit(0.5, "lines"),
        max.overlaps = Inf
      ) +
      ggtitle(paste0("Volcano Plot - ", biotype, ": ", numerator, " vs. ", denominator)) +
      xlab("log2 Fold Change") + ylab("-log10 p-value") +
      theme_minimal() +
      geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
      geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
      xlim(volcano_xlim) +
      ylim(NA, volcano_ylim_max)
    
    ggsave(file.path(results_dir, paste0("06_", biotype, "_", comp_name, "_volcano_plot.png")),
           volcano_plot, width = 8, height = 7)
  }
}

# --- Run for mRNA and lncRNA ---
run_deseq_analysis(mRNA_obj, vsd_corrected_mRNA, "mRNA")
run_deseq_analysis(lncRNA_obj, vsd_corrected_lncRNA, "lncRNA")

message("\n--- Differential expression analysis completed. Check 'results' folder for CSV and volcano plots ---")