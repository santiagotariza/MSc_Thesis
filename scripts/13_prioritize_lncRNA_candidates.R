# ============================================================
# 13_prioritize_lncRNA_candidates.R
# Filtering and prioritization of lncRNA candidates
#
# Input:
#   - 12_integrated_ORF_pathways_summary.csv
#     (integrated ORFs and pathway enrichment results)
#   - 09_lncRNA_*_coexpression_ORA_<lncRNA>.RData
#     (GO enrichment results for candidate lncRNAs, required for plots)
#
# Output:
#   - 13_top_lncRNA_candidates.csv
#     (prioritized candidate list based on filtering criteria)
#   - 13_final_candidates_dotplot_<lncRNA>.png
#     (visualizations of pathway enrichment for prioritized candidates)
#
# Description:
#   This script applies filtering criteria to prioritize lncRNAs with
#   significant ORFs and biologically relevant pathways. It also generates
#   visualization plots (dotplots) for the top candidates.
# ============================================================

# --- LOAD LIBRARIES ---
library(tidyverse)
library(stringr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# --- DIRECTORIES ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}

results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- LOAD MASTER DATA ---
integrated_data_path <- file.path(results_dir, "12_integrated_ORF_pathways_summary.csv")

if (!file.exists(integrated_data_path)) {
  stop("Integrated data file not found (12_integrated_ORF_pathways_summary.csv). Please run script 12 first.")
}

integrated_df <- read_csv(integrated_data_path, show_col_types = FALSE)

colnames(integrated_df)
class(integrated_df)
str(integrated_df)



# --- PRIORITIZATION CRITERIA ---
message("Applying filtering and prioritization criteria...")

relevance_terms <- c(
  "cancer", "tumor", "apoptosis", "proliferat",
  "metastasis", "migration", "inflammatory", "immune",
  "response to stimulus", "cell death"
)
relevance_pattern <- paste(relevance_terms, collapse = "|")


orf_counts <- dplyr::count(integrated_df, lncRNA_id, name = "n_ORFs")



# Add counts to the original dataframe
integrated_df <- left_join(integrated_df, orf_counts, by = "lncRNA_id")

# Filter
top_candidates_df <- integrated_df %>%
  filter(p.adjust < 0.05) %>%
  # From here you can filter how many ORFs should have the lncRNA
  #filter(n_ORFs <= 100) %>%
  arrange(p.adjust, desc(Count))



# --- SAVE PRIORITIZED LIST ---
output_path <- file.path(results_dir, "13_top_lncRNA_candidates.csv")
write_csv(top_candidates_df, output_path)

message(paste("\nPrioritization completed. Candidate list saved at:", output_path))
message(paste("Number of final candidates:", nrow(top_candidates_df)))

# --- GENERATE VISUALIZATIONS ---
message("\nGenerating enrichment plots for prioritized candidates...")


generate_and_save_plots <- function(lncRNA_id, comp_name, results_dir) {
  ora_results_path <- file.path(
    results_dir,
    paste0("09_lncRNA_", comp_name, "_coexpression_ORA_", lncRNA_id, ".RData")
  )
  
  if (file.exists(ora_results_path)) {
    load(ora_results_path)
    
    if (!is.null(go_enrich_lncRNA) && nrow(go_enrich_lncRNA) > 0) {
      
      #go_enrich_lncRNA$Description <- str_wrap(go_enrich_lncRNA$Description, width = 45)
      # Dotplot
      dotplot_go <- dotplot(go_enrich_lncRNA, showCategory = 10) +
        ggtitle(paste("Co-expression ORA:", lncRNA_id, "-", comp_name))
      ggsave(
        file.path(results_dir, paste0("13_final_candidates_dotplot_", lncRNA_id, ".png")),
        dotplot_go, width = 10, height = 10
      )
      
      message(paste("  - Plots generated for", lncRNA_id))
    } else {
      message(paste("  - No enrichment results for", lncRNA_id, ". Skipping plots."))
    }
  } else {
    warning(paste("  - ORA results file not found for", lncRNA_id, ". Skipping plots."))
  }
}

unique_final_candidates <- top_candidates_df %>%
  distinct(lncRNA_id, comparison)

if (nrow(unique_final_candidates) > 0) {
  for (i in 1:nrow(unique_final_candidates)) {
    lncRNA_id <- unique_final_candidates$lncRNA_id[i]
    comp_name <- unique_final_candidates$comparison[i]
    
    generate_and_save_plots(lncRNA_id, comp_name, results_dir)
  }
} else {
  message("No final candidates available for plotting.")
}

message("\nPlot generation completed.")