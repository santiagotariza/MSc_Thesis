###############################################################################
# Script: 03_normalize_transform.R
# Purpose: Create DESeq2 objects for mRNA and lncRNA, filter low-count genes, 
#          normalize counts, and apply variance stabilizing transformation (VST).
#
# Dependencies:
#   - DESeq2
#   - rstudioapi
#
# Inputs:
#   - results/01_filtered_data.RData (sample_info, counts_mRNA, counts_lncRNA)
#
# Outputs:
#   - results/03_dds_vsd_objects.RData (mRNA_obj, lncRNA_obj)
###############################################################################

# Load required libraries
library(DESeq2)
library(rstudioapi)

# Define directories
script_dir  <- dirname(rstudioapi::getSourceEditorContext()$path)
results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Load data from previous step
# -------------------------------------------------------------------------
load(file.path(results_dir, "01_filtered_data.RData"))

# -------------------------------------------------------------------------
# Function to create DESeq2 dataset and VST
# -------------------------------------------------------------------------
create_dds <- function(counts, sample_info) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts, 
    colData   = sample_info, 
    design    = ~ condition
  )
  
  # Diagnostics: dimensions before filtering
  message("Dimensions before filtering by rowSums: ", paste(dim(dds), collapse = " x "))
  
  # Filter out low-count genes
  dds <- dds[rowSums(counts(dds)) > 100, ]
  
  # Diagnostics: dimensions after filtering
  message("Dimensions after filtering by rowSums: ", paste(dim(dds), collapse = " x "))
  
  # Check if dataset is empty after filtering
  if (nrow(dds) == 0) {
    warning("DESeqDataSet is empty after filtering by rowSums. Cannot proceed with DESeq and VST.")
    return(NULL)
  }
  
  # Run DESeq and variance stabilizing transformation
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)
  
  return(list(dds = dds, vsd = vsd))
}

# -------------------------------------------------------------------------
# Process mRNA
# -------------------------------------------------------------------------
message("Processing mRNA:")
mRNA_obj <- create_dds(counts_mRNA, sample_info)

# -------------------------------------------------------------------------
# Process lncRNA
# -------------------------------------------------------------------------
message("Processing lncRNA:")
lncRNA_obj <- create_dds(counts_lncRNA, sample_info)

# -------------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------------
if (!is.null(mRNA_obj) && !is.null(lncRNA_obj)) {
  save(mRNA_obj, lncRNA_obj, file = file.path(results_dir, "03_dds_vsd_objects.RData"))
  message("Saved DESeq2 and VST objects to 03_dds_vsd_objects.RData")
} else {
  warning("Could not create DESeq2/VST objects for mRNA or lncRNA. Skipping save step.")
}