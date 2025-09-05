# ============================================================
# 14_render_rna_seq_report.R
# Script to render the final RNA-seq report as HTML
#
# Input:
#   - rna_seq_report.Rmd
#       (R Markdown file with all analysis results and plots)
#
# Output:
#   - results/rna_seq_report.html
#       (HTML report summarizing the RNA-seq analysis)
#
# Description:
#   This script renders the R Markdown report into an HTML document.
#   It ensures that relative paths work correctly by setting the working
#   directory to the script location, and saves the report in the results folder.
# ============================================================

# --- LOAD LIBRARIES ---
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  install.packages("rmarkdown")
}
library(rmarkdown)

# --- DIRECTORIES ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}

results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

rmd_path <- file.path(script_dir, "rna_seq_report.Rmd")
output_file <- file.path(results_dir, "rna_seq_report.html")

# --- RENDER REPORT ---
message("Rendering the RNA-seq report...")

# Set working directory to script location to ensure relative paths in Rmd work
setwd(script_dir)

render(
  input = basename(rmd_path),          # Use only the filename; wd is now correct
  output_format = "html_document",
  output_dir = results_dir,
  output_file = "rna_seq_report.html",
  quiet = TRUE
)

message(paste("Report rendered successfully and saved at:", output_file))