# ============================================================
# 10_consolidate_coexp_results.R
# Consolidation of lncRNA co-expression ORA results
#
# Input:
#   - 09_lncRNA_*_coexpression_ORA_results_<lncRNA>.csv
#       (ORA results from script 09, one file per lncRNA per comparison)
#
# Output:
#   - 10_consolidated_coexpression_ORA_results.csv
#
# Description:
#   This script consolidates individual lncRNA co-expression ORA
#   results into a single CSV file. It extracts comparison names and
#   lncRNA IDs from file names, merges all data, and saves the
#   consolidated results.
# ============================================================

# --- LOAD LIBRARIES ---
library(tidyverse)
library(stringr)

# --- DIRECTORIES ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}

results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- GET LIST OF FILES ---
file_list <- list.files(
  results_dir,
  pattern = "^09_lncRNA_.*_coexpression_ORA_results_.*\\.csv$",
  full.names = TRUE
)

if (length(file_list) == 0) {
  message("No ORA co-expression result files found. Run script 09 first.")
} else {
  message(paste("Found", length(file_list), "files to consolidate."))
  
  # --- CONSOLIDATE RESULTS ---
  consolidated_df <- tibble()
  
  for (file_path in file_list) {
    file_name <- basename(file_path)
    # Example: 09_lncRNA_G3_HMGB1_KO_vs_G1_NTC_coexpression_ORA_results_ENSG00000287472.csv
    
    # Extract lncRNA ID
    parts <- str_split(file_name, "_results_")[[1]]
    lncRNA_id <- str_replace(parts[2], "\\.csv$", "")
    
    # Extract comparison name
    comp_part <- str_replace(parts[1], "^09_lncRNA_", "")
    comp_name <- str_replace(comp_part, "_coexpression_ORA", "")
    
    # Read CSV
    temp_df <- read_csv(file_path, show_col_types = FALSE) %>%
      mutate(
        lncRNA_id = lncRNA_id,
        comparison = comp_name
      )
    
    # Append to consolidated dataframe
    consolidated_df <- bind_rows(consolidated_df, temp_df)
  }
  
  # Reorder columns
  consolidated_df <- consolidated_df %>%
    relocate(lncRNA_id, comparison, .before = 1)
  
  # --- SAVE CONSOLIDATED FILE ---
  output_path <- file.path(results_dir, "10_consolidated_coexpression_ORA_results.csv")
  write_csv(consolidated_df, output_path)
  
  message(paste("\nConsolidation completed. Results saved at:", output_path))
}