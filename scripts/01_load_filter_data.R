###############################################################################
# Script: 01_load_filter_data.R
# Purpose: Load raw gene counts, clean sample names, assign conditions, 
#          extract gene biotypes from GTF, filter by type (mRNA and lncRNA),
#          and save processed objects for downstream analysis.
#
# Dependencies:
#   - tidyverse
#   - rtracklayer
#   - DESeq2
#   - rstudioapi
#
# Inputs:
#   - counts/gene_counts.txt
#   - genome/gencode.v48.primary_assembly.annotation.gtf
#
# Outputs:
#   - results/01_filtered_data.RData
###############################################################################

# Load required libraries
library(tidyverse)
library(rtracklayer)
library(DESeq2)
library(rstudioapi)

# Define file paths
script_dir   <- dirname(rstudioapi::getSourceEditorContext()$path)
counts_file  <- file.path(script_dir, "../counts/gene_counts.txt")
gtf_file     <- file.path(script_dir, "../genome/gencode.v48.primary_assembly.annotation.gtf")
results_dir  <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load counts and metadata
# -------------------------------------------------------------------------
message("Loading counts data...")
counts_data <- read.delim(
  counts_file, 
  comment.char = "#", 
  row.names = 1, 
  check.names = FALSE
)

# Extract raw column names before R modifies them
raw_names   <- colnames(counts_data)[6:ncol(counts_data)]
clean_names <- sub("\\.Aligned.*", "", basename(raw_names))

# Build count matrix with cleaned sample names
counts_matrix <- counts_data[, 6:ncol(counts_data)]
colnames(counts_matrix) <- clean_names
rownames(counts_matrix) <- sub("\\..*", "", rownames(counts_matrix))

# -------------------------------------------------------------------------
# Build sample metadata
# -------------------------------------------------------------------------
message("Building sample metadata...")
sample_info <- data.frame(
  sample = colnames(counts_matrix),
  stringsAsFactors = FALSE
)

# Assign conditions based on sample name prefixes
sample_info$condition <- factor(
  case_when(
    grepl("^NTC", sample_info$sample)   ~ "G1_NTC",
    grepl("^WTC", sample_info$sample)   ~ "G2_WTC",
    grepl("^H1KOC", sample_info$sample) ~ "G3_HMGB1_KO",
    grepl("^H2KOC", sample_info$sample) ~ "G4_HMGB2_KO",
    TRUE ~ "Unknown"
  ),
  levels = c("G1_NTC", "G2_WTC", "G3_HMGB1_KO", "G4_HMGB2_KO")
)

rownames(sample_info) <- sample_info$sample
sample_info <- sample_info[match(colnames(counts_matrix), rownames(sample_info)), ]

# -------------------------------------------------------------------------
# Load GTF file and extract gene biotypes
# -------------------------------------------------------------------------
message("Importing GTF file...")
gtf <- import(gtf_file)
gtf_genes <- gtf[gtf$type == "gene"]

if ("gene_type" %in% names(mcols(gtf_genes))) {
  gene_biotypes <- mcols(gtf_genes)[, c("gene_id", "gene_type")]
} else if ("gene_biotype" %in% names(mcols(gtf_genes))) {
  gene_biotypes <- mcols(gtf_genes)[, c("gene_id", "gene_biotype")]
  names(gene_biotypes)[2] <- "gene_type"
} else {
  stop("No gene_type or gene_biotype found in GTF.")
}

# Clean and deduplicate gene biotypes
gene_biotypes$gene_id <- sub("\\..*", "", gene_biotypes$gene_id)
gene_biotypes <- unique(as.data.frame(gene_biotypes))

# -------------------------------------------------------------------------
# Filter genes by biotype
# -------------------------------------------------------------------------
message("Filtering genes by biotype...")
mRNA_genes  <- gene_biotypes$gene_id[gene_biotypes$gene_type %in% c("protein_coding")]
lncRNA_genes <- gene_biotypes$gene_id[gene_biotypes$gene_type %in% c(
  "lncRNA", "antisense", "sense_intronic", "sense_overlapping", 
  "3prime_overlapping_ncrna"
)]

counts_mRNA   <- counts_matrix[rownames(counts_matrix) %in% mRNA_genes, ]
counts_lncRNA <- counts_matrix[rownames(counts_matrix) %in% lncRNA_genes, ]

# -------------------------------------------------------------------------
# Save processed objects
# -------------------------------------------------------------------------
save(
  sample_info, 
  counts_mRNA, 
  counts_lncRNA, 
  gene_biotypes, 
  file = file.path(results_dir, "01_filtered_data.RData")
)

message("Step 01 completed. Filtered data saved to 01_filtered_data.RData")