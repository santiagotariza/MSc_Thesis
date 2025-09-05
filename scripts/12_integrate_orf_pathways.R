# ============================================================
# 12_integrate_orf_pathways.R
# Integration of predicted ORFs with co-expression pathway results
#
# Input:
#   - 11_predicted_lncRNA_ORFs.csv
#       (predicted ORFs from script 11)
#   - 10_consolidated_coexpression_ORA_results.csv
#       (consolidated ORA results from script 10)
#
# Output:
#   - 12_integrated_ORF_pathways_summary.csv
#       (summary table linking ORFs with co-expression pathways)
#
# Description:
#   This script integrates the predicted ORFs from lncRNAs with their
#   corresponding co-expression pathway enrichment results, generating
#   a combined summary for downstream interpretation.
# ============================================================

# --- LOAD LIBRARIES ---
library(tidyverse)
library(dplyr)

# --- DIRECTORIES ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}

results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- LOAD DATA ---
message("Loading ORF and co-expression pathway results...")

orf_results_path <- file.path(results_dir, "11_predicted_lncRNA_ORFs.csv")
coexp_results_path <- file.path(results_dir, "10_consolidated_coexpression_ORA_results.csv")

if (!file.exists(orf_results_path) || !file.exists(coexp_results_path)) {
  stop("Required input files not found. Please run scripts 10 and 11 first.")
}

orf_df <- read_csv(orf_results_path, show_col_types = FALSE)
coexp_df <- read_csv(coexp_results_path, show_col_types = FALSE)

# --- MERGE DATA ---
message("Merging ORF and pathway results...")

coexp_slim <- coexp_df %>%
  dplyr::select(lncRNA_id, comparison, Description, p.adjust, Count) %>%
  distinct(lncRNA_id, Description, .keep_all = TRUE)  # remove duplicate entries

final_integrated_df <- left_join(
  orf_df,
  coexp_slim,
  by = "lncRNA_id",
  relationship = "many-to-many"
)

final_integrated_df <- final_integrated_df %>%
  relocate(comparison, .after = lncRNA_id) %>%
  arrange(comparison, lncRNA_id, start)

# --- SAVE OUTPUT ---
output_path <- file.path(results_dir, "12_integrated_ORF_pathways_summary.csv")
write_csv(final_integrated_df, output_path)

message(paste("\nIntegration completed. Final summary saved at:", output_path))