#!/usr/bin/env Rscript

# =========================================================
# 00_run_all_analysis.R - Executes the complete RNA-seq pipeline
# =========================================================

# --- DIRECTORIES ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}

setwd(script_dir)  # <- Añade esta línea

# 1. Guarda el tiempo de inicio
tiempo_inicio <- proc.time()




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
cat("STARTING COMPLETE RNA-seq PIPELINE\n")
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
cat("RNA-seq PIPELINE COMPLETED SUCCESSFULLY\n")
cat("Results available in: results/\n")
cat("Fastp reports: fastp_reports/\n")
cat("Final HTML report: results/rna_seq_report.html\n")
cat("====================================================\n")

# 3. Guarda el tiempo de finalización
tiempo_final <- proc.time()

# 4. Calcula y muestra la diferencia
tiempo_ejecucion <- tiempo_final - tiempo_inicio
print(tiempo_ejecucion)

# install.packages("renv")
# 
# renv::init()  # solo la primera vez
# 
# 
# renv::snapshot()
