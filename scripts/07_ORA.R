##############################################
# Script: 07_ORA.R
# Purpose: Perform Over-Representation Analysis (ORA) for DE genes
# Input:   DE results (from script 06)
# Output:  ORA results (RData, CSV), dotplots (PNG)
##############################################

# --- Libraries ---
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(enrichplot)
library(ReactomePA)

# --- Parameters ---
log2FC_threshold <- 1
padj_threshold <- 0.01

# --- Define directories ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}
results_dir <- normalizePath(file.path(script_dir, "..", "results"), mustWork = FALSE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper functions ---
gene_id_to_entrez <- function(gene_id_list) {
  bitr(gene_id_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
}

# --- Define comparisons of interest ---
comparisons <- c(
  "G3_HMGB1_KO_vs_G1_NTC",
  "G4_HMGB2_KO_vs_G1_NTC",
  "G2_WTC_vs_G1_NTC",
  "G3_HMGB1_KO_vs_G4_HMGB2_KO",
  "G3_HMGB1_KO_vs_G2_WTC",
  "G4_HMGB2_KO_vs_G2_WTC"
)

# --- Run ORA for mRNA genes ---
message("\n--- Running ORA for mRNA genes ---")

for (comp_name in comparisons) {
  message(paste("  - Processing ORA for:", comp_name))
  
  de_results_path <- file.path(results_dir, paste0("06_mRNA_", comp_name, "_DE_results.csv"))
  if (!file.exists(de_results_path)) {
    warning(paste("File not found:", de_results_path, ". Skipping ORA for this comparison."))
    next
  }
  
  de_results <- read.csv(de_results_path, stringsAsFactors = FALSE)
  significant_genes <- de_results$gene_id
  entrez_ids <- gene_id_to_entrez(significant_genes)
  
  if (nrow(entrez_ids) > 0) {
    # 1. ORA with Gene Ontology (GO)
    go_enrich <- clusterProfiler::enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    
    # 2. ORA with Reactome pathways
    reactome_enrich <- ReactomePA::enrichPathway(
      gene = entrez_ids$ENTREZID,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    
    # Save ORA objects (RData)
    ora_output_path <- file.path(results_dir, paste0("07_mRNA_ORA_", comp_name, ".RData"))
    save(go_enrich, reactome_enrich, file = ora_output_path)
    message(paste("    - ORA completed. RData results saved to:", ora_output_path))
    
    # Save and plot GO results
    if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
      write.csv(as.data.frame(go_enrich), 
                file = file.path(results_dir, paste0("07_mRNA_ORA_GO_", comp_name, ".csv")), 
                row.names = FALSE)
      
      dotplot_go <- dotplot(go_enrich, showCategory = 15) +
        ggtitle(paste("ORA GO Dotplot:", comp_name))
      ggsave(file.path(results_dir, paste0("07_mRNA_ORA_GO_dotplot_", comp_name, ".png")),
             dotplot_go, width = 10, height = 8)
    } else {
      message("    - No significant GO results found. Skipping.")
    }
    
    # Save and plot Reactome results
    if (!is.null(reactome_enrich) && nrow(as.data.frame(reactome_enrich)) > 0) {
      write.csv(as.data.frame(reactome_enrich), 
                file = file.path(results_dir, paste0("07_mRNA_ORA_Reactome_", comp_name, ".csv")), 
                row.names = FALSE)
      
      dotplot_reactome <- dotplot(reactome_enrich, showCategory = 15) +
        ggtitle(paste("ORA Reactome Dotplot:", comp_name))
      ggsave(file.path(results_dir, paste0("07_mRNA_ORA_Reactome_dotplot_", comp_name, ".png")),
             dotplot_reactome, width = 10, height = 8)
    } else {
      message("    - No significant Reactome results found. Skipping.")
    }
    
  } else {
    message("    - Could not map significant genes to EntrezID. Skipping ORA.")
  }
}

message("\n--- ORA analysis completed ---")