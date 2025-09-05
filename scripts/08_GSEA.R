##############################################
# Script: 08_GSEA.R
# Purpose: Perform Gene Set Enrichment Analysis (GSEA) for DE genes
# Input:   DE results (from script 06) and dds objects
# Output:  GSEA results (RData, CSV), dotplots and enrichment plots (PNG)
##############################################

# --- Libraries ---
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(enrichplot)
library(ReactomePA)
library(DESeq2)

# --- Parameters ---
padj_threshold <- 0.01

# --- Define directories ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}
results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- Define comparisons of interest ---
comparisons <- c(
  "G3_HMGB1_KO_vs_G1_NTC",
  "G4_HMGB2_KO_vs_G1_NTC",
  "G2_WTC_vs_G1_NTC",
  "G3_HMGB1_KO_vs_G4_HMGB2_KO",
  "G3_HMGB1_KO_vs_G2_WTC",
  "G4_HMGB2_KO_vs_G2_WTC"
)

# --- Load expression data ---
load(file.path(results_dir, "03_dds_vsd_objects.RData"))

# --- Run GSEA for mRNA genes ---
message("\n--- Running GSEA for mRNA genes ---")

for (comp_name in comparisons) {
  message(paste("  - Processing GSEA for:", comp_name))
  
  de_results_path <- file.path(results_dir, paste0("06_mRNA_", comp_name, "_DE_results.csv"))
  if (!file.exists(de_results_path)) {
    warning(paste("File not found:", de_results_path, ". Skipping GSEA for this comparison."))
    next
  }
  
  dds <- mRNA_obj$dds
  
  # --- Shrink log2FC (coef vs contrast) ---
  comp_parts <- str_split(comp_name, "_vs_")[[1]]
  group_1 <- comp_parts[1]
  group_2 <- comp_parts[2]
  
  res_shrunk <- if (str_detect(comp_name, "G1_NTC$")) {
    coef_name <- paste0("condition_", group_1, "_vs_G1_NTC")
    
    tryCatch({
      lfcShrink(dds, coef = coef_name, type = "apeglm")
    }, error = function(e) {
      message(paste("    Error with apeglm:", e$message, "Trying with ashr."))
      lfcShrink(dds, coef = coef_name, type = "ashr")
    })
    
  } else {
    res_raw <- results(dds, contrast = c("condition", group_1, group_2))
    
    tryCatch({
      lfcShrink(dds, contrast = c("condition", group_1, group_2), type = "ashr")
    }, error = function(e) {
      message(paste("    Error with lfcShrink:", e$message, "Using raw results without shrinkage."))
      as.data.frame(res_raw)
    })
  }
  
  # --- Prepare ranked gene list for GSEA ---
  gene_list_df <- data.frame(ENSEMBL = rownames(res_shrunk), log2FoldChange = res_shrunk$log2FoldChange)
  
  gene_list_df <- gene_list_df %>%
    left_join(bitr(unique(gene_list_df$ENSEMBL), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db), by = "ENSEMBL") %>%
    filter(!is.na(ENTREZID)) %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1) %>%
    ungroup() %>%
    dplyr::select(ENTREZID, log2FoldChange)
  
  final_gene_list <- gene_list_df$log2FoldChange
  names(final_gene_list) <- gene_list_df$ENTREZID
  final_gene_list <- sort(final_gene_list, decreasing = TRUE)
  
  if (length(final_gene_list) > 0) {
    # 1. GSEA with Gene Ontology (GO)
    gsea_go <- clusterProfiler::gseGO(
      geneList = final_gene_list,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pvalueCutoff = 0.05,
      verbose = FALSE
    )
    
    # 2. GSEA with Reactome pathways
    gsea_reactome <- ReactomePA::gsePathway(
      geneList = final_gene_list,
      pvalueCutoff = 0.05,
      verbose = FALSE
    )
    
    # Save GSEA results (RData)
    gsea_output_path <- file.path(results_dir, paste0("08_mRNA_GSEA_", comp_name, ".RData"))
    save(gsea_go, gsea_reactome, file = gsea_output_path)
    message(paste("    - GSEA completed. RData results saved to:", gsea_output_path))
    
    # --- Save and plot GO results ---
    if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
      write.csv(as.data.frame(gsea_go), 
                file = file.path(results_dir, paste0("08_mRNA_GSEA_GO_", comp_name, ".csv")), 
                row.names = FALSE)
      
      dotplot_gsea_go <- dotplot(gsea_go, showCategory = 15) +
        ggtitle(paste("GSEA GO Dotplot:", comp_name))
      ggsave(file.path(results_dir, paste0("08_mRNA_GSEA_GO_dotplot_", comp_name, ".png")),
             dotplot_gsea_go, width = 10, height = 8)
      
      gsea_plot <- gseaplot2(gsea_go, geneSetID = 1:min(5, nrow(gsea_go)),
                             title = paste("GSEA - Top Pathways for mRNA:", comp_name))
      ggsave(file.path(results_dir, paste0("08_mRNA_GSEA_gseaplot2_", comp_name, ".png")),
             gsea_plot, width = 10, height = 8)
    } else {
      message("    - No significant GO results found. Skipping.")
    }
    
    # --- Save and plot Reactome results ---
    if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) > 0) {
      write.csv(as.data.frame(gsea_reactome), 
                file = file.path(results_dir, paste0("08_mRNA_GSEA_Reactome_", comp_name, ".csv")), 
                row.names = FALSE)
      
      dotplot_reactome_gsea <- dotplot(gsea_reactome, showCategory = 15) +
        ggtitle(paste("GSEA Reactome Dotplot:", comp_name))
      ggsave(file.path(results_dir, paste0("08_mRNA_GSEA_Reactome_dotplot_", comp_name, ".png")),
             dotplot_reactome_gsea, width = 10, height = 8)
    } else {
      message("    - No significant Reactome results found. Skipping.")
    }
  } else {
    message("    - Could not map genes for GSEA. Skipping analysis.")
  }
}

message("\n--- GSEA analysis completed ---")