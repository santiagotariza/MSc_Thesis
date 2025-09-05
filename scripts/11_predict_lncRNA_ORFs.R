# ============================================================
# 11_predict_lncRNA_ORFs.R
# Prediction of ORFs in differentially expressed lncRNAs
#
# Input:
#   - 06_lncRNA_<comparison>_DE_results.csv
#       (DE results for lncRNAs from script 06)
#   - Ensembl (via biomaRt) → transcript sequences (cDNA) of lncRNAs
#
# Output:
#   - 11_predicted_lncRNA_ORFs.csv
#       (table with ORFs predicted from lncRNAs)
#   - 11_predicted_lncRNA_ORF_sequences.fasta
#       (FASTA file with ORF nucleotide sequences)
#
# Description:
#   This script gathers all significantly differentially expressed
#   lncRNAs, retrieves their transcript sequences from Ensembl,
#   predicts open reading frames (ORFs) using ORFik, and outputs
#   a summary table and FASTA sequences for ORFs above a minimum
#   length threshold.
# ============================================================

# Instalar BiocManager si no está instalado
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar paquetes sin interacción
BiocManager::install("biomaRt", ask = FALSE)
BiocManager::install("ORFik", ask = FALSE)


# --- LOAD LIBRARIES ---
library(tidyverse)
library(biomaRt)
library(Biostrings)
library(stringr)
library(ORFik)
library(dplyr)

# --- PARAMETERS ---
padj_threshold <- 0.01
min_orf_length_aa <- 50  # minimum ORF length in amino acids

# --- DIRECTORIES ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}

results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- GET CANDIDATE lncRNAs ---
lncRNA_comparisons <- c(
  "G3_HMGB1_KO_vs_G1_NTC",
  "G4_HMGB2_KO_vs_G1_NTC",
  "G2_WTC_vs_G1_NTC",
  "G3_HMGB1_KO_vs_G4_HMGB2_KO",
  "G3_HMGB1_KO_vs_G2_WTC",
  "G4_HMGB2_KO_vs_G2_WTC"
)

all_de_lncRNAs <- tibble()

for (comp_name in lncRNA_comparisons) {
  de_path <- file.path(results_dir, paste0("06_lncRNA_", comp_name, "_DE_results.csv"))
  
  if (file.exists(de_path)) {
    de_results <- read_csv(
      de_path,
      col_types = cols(
        gene_id = col_character(),
        baseMean = col_double(),
        log2FoldChange = col_double(),
        lfcSE = col_double(),
        stat = col_double(),
        pvalue = col_double(),
        padj = col_double()
      )
    ) %>%
      filter(padj < padj_threshold) %>%
      mutate(comparison = comp_name)
    
    all_de_lncRNAs <- bind_rows(all_de_lncRNAs, de_results)
  } else {
    warning(paste("File not found:", de_path))
  }
}

unique_lncRNA_ids <- unique(all_de_lncRNAs$gene_id)

if (length(unique_lncRNA_ids) == 0) {
  stop("No significant lncRNAs found. Stopping process.")
}
message(paste("Found", length(unique_lncRNA_ids), "significant lncRNAs for ORF analysis."))

# --- RETRIEVE SEQUENCES FROM ENSEMBL ---
message("Connecting to Ensembl to retrieve lncRNA sequences...")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

sequences <- getSequence(
  id = unique_lncRNA_ids,
  type = "ensembl_gene_id",
  seqType = "cdna",
  mart = ensembl
)

# Filter out missing sequences
sequences <- sequences %>% drop_na(cdna)
unique_lncRNA_ids <- unique(sequences$ensembl_gene_id)
message(paste("Retrieved sequences for", length(unique_lncRNA_ids), "lncRNAs."))

# --- ORF PREDICTION ---
message("Predicting ORFs in lncRNA sequences...")

lncRNA_sequences <- DNAStringSet(sequences$cdna)
names(lncRNA_sequences) <- sequences$ensembl_gene_id

all_orfs <- findORFs(lncRNA_sequences, startCodon = "ATG")
all_orfs_df <- as.data.frame(all_orfs)

# --- PROCESS RESULTS ---
message("Processing ORF prediction results...")

predicted_orfs_df <- all_orfs_df %>%
  as_tibble() %>%
  mutate(
    length_nt = width,
    length_aa = width / 3,
    lncRNA_id = names(lncRNA_sequences)[group]
  ) %>%
  filter(length_aa >= min_orf_length_aa) %>%
  arrange(lncRNA_id, desc(length_aa))

message(paste(
  "Predicted", nrow(predicted_orfs_df),
  "ORFs with length >=", min_orf_length_aa, "aa."
))

# --- SAVE RESULTS ---
# 1. Save summary table
output_csv_path <- file.path(results_dir, "11_predicted_lncRNA_ORFs.csv")
write_csv(predicted_orfs_df, output_csv_path)
message(paste("ORF summary table saved at:", output_csv_path))

# 2. Save ORF sequences in FASTA format
if (nrow(predicted_orfs_df) > 0) {
  orf_sequences <- DNAStringSet()
  
  for (i in seq_len(nrow(predicted_orfs_df))) {
    row <- predicted_orfs_df[i, ]
    seq <- lncRNA_sequences[row$lncRNA_id]
    
    if (length(seq) > 0 && row$start >= 1 && row$end <= width(seq)) {
      orf_seq <- subseq(seq, start = row$start, end = row$end)
      
      names(orf_seq) <- paste0(
        row$lncRNA_id,
        "|start:", row$start,
        "|end:", row$end,
        "|len_nt:", row$length_nt,
        "|len_aa:", row$length_aa
      )
      
      orf_sequences <- c(orf_sequences, orf_seq)
    } else {
      warning(paste("Invalid coordinates for ORF in", row$lncRNA_id))
    }
  }
  
  output_fasta_path <- file.path(results_dir, "11_predicted_lncRNA_ORF_sequences.fasta")
  writeXStringSet(orf_sequences, output_fasta_path)
  message(paste("ORF sequences saved in FASTA format at:", output_fasta_path))
}

message("\nORF prediction analysis completed.")