# ============================================================
# 09_coexp_lncRNA.R
# lncRNA Co-expression Analysis
#
# Input:
#   - 03_dds_vsd_objects.RData : contains mRNA_obj and lncRNA_obj (with vsd)
#   - 06_lncRNA_*_DE_results.csv : differential expression results for lncRNAs
#
# Output:
#   - 09_lncRNA_*_coexpression_ORA_<lncRNA>.RData
#   - 09_lncRNA_*_coexpression_ORA_results_<lncRNA>.csv
#
# Description:
#   For each comparison, the top significant lncRNAs are selected.
#   Co-expression with mRNAs is computed using Pearson correlation.
#   The top correlated mRNAs are mapped to Entrez IDs and subjected
#   to Over-Representation Analysis (ORA) using Gene Ontology (GO).
# ============================================================

# --- LOAD LIBRARIES ---
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(enrichplot)
library(ggplot2)


# --- PARAMETERS ---
padj_threshold <- 0.05
log2fc_threshold <- 1
baseMean_threshold <- 10
top_lncRNAs_to_analyze <- 5 


# --- DIRECTORIES ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}

results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- HELPER FUNCTIONS ---
gene_id_to_entrez <- function(gene_id_list) {
  bitr(
    gene_id_list,
    fromType = "ENSEMBL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
  )
}

# --- LOAD DATA ---
load(file.path(results_dir, "03_dds_vsd_objects.RData"))
vsd_mRNA   <- mRNA_obj$vsd
vsd_lncRNA <- lncRNA_obj$vsd

# --- COMPARISONS ---
lncRNA_comparisons <- c(
  "G3_HMGB1_KO_vs_G1_NTC",
  "G4_HMGB2_KO_vs_G1_NTC",
  "G2_WTC_vs_G1_NTC"
)

# --- RUNNING CO-EXPRESSION ANALYSIS ---
message("\n--- Running co-expression analysis for lncRNAs ---")

for (comp_name in lncRNA_comparisons) {
  message(paste("\n--- Processing comparison:", comp_name, "---"))
  
  lncRNA_de_results_path <- file.path(
    results_dir,
    paste0("06_lncRNA_", comp_name, "_DE_results.csv")
  )
  
  if (file.exists(lncRNA_de_results_path)) {
    de_results <- read.csv(lncRNA_de_results_path, stringsAsFactors = FALSE)
    
    if (nrow(de_results) > 0) {
      # Select top lncRNAs for this comparison
      top_lncRNAs <- de_results %>%
        filter(
          padj < padj_threshold,
          abs(log2FoldChange) > log2fc_threshold,
          baseMean > baseMean_threshold
        ) %>%
        mutate(relevance_score = -log10(padj) * abs(log2FoldChange)) %>%
        arrange(desc(relevance_score)) %>%
        head(top_lncRNAs_to_analyze)
      
      message(paste("    -", nrow(top_lncRNAs), "selected lncRNAs after applying filters"))
      
      
      
      if (nrow(top_lncRNAs) > 0) {
        for (i in 1:nrow(top_lncRNAs)) {
          lncRNA_id <- top_lncRNAs$gene_id[i]
          message(paste("  - Analyzing co-expression for lncRNA:", lncRNA_id, 
                        "(padj =", round(top_lncRNAs$padj[i], 4), ")"))
          
          # Pearson correlation between lncRNA and all mRNAs
          cor_values <- cor(
            t(assay(vsd_lncRNA)[lncRNA_id, , drop = FALSE]),
            t(assay(vsd_mRNA)),
            method = "pearson"
          )
          cor_df <- data.frame(
            gene_id     = colnames(cor_values),
            correlation = as.vector(cor_values)
          )
          
          # Select most correlated mRNAs
          top_correlated_mRNAs <- cor_df %>%
            arrange(desc(abs(correlation))) %>%
            head(200)
          
          # Map to Entrez IDs and run ORA
          entrez_ids_mRNAs <- gene_id_to_entrez(top_correlated_mRNAs$gene_id)
          
          if (nrow(entrez_ids_mRNAs) > 0) {
            go_enrich_lncRNA <- clusterProfiler::enrichGO(
              gene          = entrez_ids_mRNAs$ENTREZID,
              OrgDb         = org.Hs.eg.db,
              ont           = "BP",
              pAdjustMethod = "BH",
              qvalueCutoff  = 0.05,
              readable      = TRUE
            )
            
            # Save RData
            coexp_output_path <- file.path(
              results_dir,
              paste0("09_lncRNA_", comp_name, "_coexpression_ORA_", lncRNA_id, ".RData")
            )
            save(go_enrich_lncRNA, file = coexp_output_path)
            message(paste("    - Co-expression ORA completed. Saved:", coexp_output_path))
            
            # Save CSV
            if (!is.null(go_enrich_lncRNA) && nrow(go_enrich_lncRNA) > 0) {
              results_df <- as.data.frame(go_enrich_lncRNA)
              csv_output_path <- file.path(
                results_dir,
                paste0("09_lncRNA_", comp_name, "_coexpression_ORA_results_", lncRNA_id, ".csv")
              )
              write.csv(results_df, file = csv_output_path, row.names = FALSE)
              message(paste("    - Co-expression ORA results saved as CSV:", csv_output_path))
            }
            
          } else {
            message("    - Could not map co-expressed mRNAs to Entrez IDs. Skipping ORA.")
          }
        }
      } else {
        message(paste("    - No significant lncRNAs found in comparison:", comp_name, ". Skipping."))
      }
    } else {
      message(paste("    - Empty results in comparison:", comp_name, ". Skipping."))
    }
  } else {
    warning(paste("    - File not found:", lncRNA_de_results_path, ". Skipping analysis."))
  }
}

message("\nCo-expression analysis completed.")