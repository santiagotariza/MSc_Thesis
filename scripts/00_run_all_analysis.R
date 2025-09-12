#!/usr/bin/env Rscript

# =========================================================
# 00_run_all_analysis.R - Executes the complete RNA-seq Analysis pipeline
# =========================================================

# Estimated processing time = 15min

# --- DIRECTORIES ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}

setwd(script_dir)


# List of scripts in execution order
scripts <- c(
  "01_load_filter_data.R",
  "02_exploratory_qc_analysis.R",
  "03_normalize_transform.R",
  "04_exploratory_plots.R",
  "05_correct_batch_effect.R",
  "06_differential_expression.R",
  "07_ORA.R",
  "08_GSEA.R",
  "09_coexp_lncRNA_groups.R",
  "10_consolidate_coexp_results.R",
  "11_predict_lncRNA_ORFs.R",
  "12_integrate_orf_pathways.R",
  "13_prioritize_lncRNA_candidates.R",
  "14_render_report.R"
)

cat("====================================================\n")
cat("STARTING RNA-seq ANALYSIS\n")
cat("====================================================\n\n")

for (script in scripts) {
  cat(paste0("Running: ", script, " ...\n"))
  
  if (!file.exists(script)) {
    stop(paste0("ERROR: Script not found: ", script))
  }
  
  tryCatch({
    source(script)
    cat(paste0("✔ Successfully executed: ", script, "\n\n"))
  }, error = function(e) {
    cat(paste0("✖ Error while running ", script, ":\n"))
    cat(e$message, "\n")
    stop("Pipeline stopped due to errors.")
  })
}

cat("====================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("Results available in: results/\n")
cat("Fastp reports: fastp_reports/\n")
cat("Final HTML report: results/rna_seq_report.html\n")
cat("====================================================\n")
