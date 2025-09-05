# 09_coexp_lncRNA.R - Script para análisis de co-expresión de lncRNA

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(enrichplot) # Se añade para los gráficos
library(ggplot2)

# --- CONFIGURACIÓN DE PARÁMETROS ---
padj_threshold <- 0.05
top_lncRNAs_to_analyze <- 5 # Define el número de lncRNAs principales a analizar para la co-expresión.

# --- Definir directorios ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}

results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- FUNCIONES AUXILIARES ---
gene_id_to_entrez <- function(gene_id_list) {
  bitr(gene_id_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
}

# Cargar los datos de lncRNA y mRNA vsd
load(file.path(results_dir, "03_dds_vsd_objects.RData"))
vsd_mRNA <- mRNA_obj$vsd
vsd_lncRNA <- lncRNA_obj$vsd

# --- REALIZANDO ANÁLISIS DE CO-EXPRESIÓN PARA lncRNA ---
message("\n--- Realizando análisis de co-expresión para lncRNAs ---")

# Obtener una lista maestra de lncRNAs significativos de todas las comparaciones
all_lncRNAs_significant <- data.frame()

# Definir las comparaciones de interés
lncRNA_comparisons <- c(
  "G3_HMGB1_KO_vs_G1_NTC",
  "G4_HMGB2_KO_vs_G1_NTC",
  "G2_WTC_vs_G1_NTC",
  "G3_HMGB1_KO_vs_G4_HMGB2_KO",
  "G3_HMGB1_KO_vs_G2_WTC",
  "G4_HMGB2_KO_vs_G2_WTC"
)
for (comp_name in lncRNA_comparisons) {
  lncRNA_de_results_path <- file.path(results_dir, paste0("06_lncRNA_", comp_name, "_DE_results.csv"))
  
  if (file.exists(lncRNA_de_results_path)) {
    de_results <- read.csv(lncRNA_de_results_path, stringsAsFactors = FALSE)
    if (nrow(de_results) > 0) {
      de_results$comparison <- comp_name
      all_lncRNAs_significant <- bind_rows(all_lncRNAs_significant, de_results)
    }
  }
}

if (nrow(all_lncRNAs_significant) > 0) {
  top_lncRNAs <- all_lncRNAs_significant %>%
    group_by(gene_id) %>%
    arrange(padj) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    arrange(padj) %>%
    head(top_lncRNAs_to_analyze)
  
  if(nrow(top_lncRNAs) > 0) {
    for (i in 1:nrow(top_lncRNAs)) {
      lncRNA_id <- top_lncRNAs$gene_id[i]
      message(paste("  - Analizando co-expresión para lncRNA:", lncRNA_id, "(padj = ", round(top_lncRNAs$padj[i], 4), " de la comparación ", top_lncRNAs$comparison[i], ")"))
      
      # Calcular correlación de Pearson entre el lncRNA y todos los mRNAs
      cor_values <- cor(t(assay(vsd_lncRNA)[lncRNA_id, , drop = FALSE]), t(assay(vsd_mRNA)), method = "pearson")
      cor_df <- data.frame(gene_id = colnames(cor_values), correlation = as.vector(cor_values))
      
      # Tomar los mRNAs más fuertemente correlacionados
      top_correlated_mRNAs <- cor_df %>%
        arrange(desc(abs(correlation))) %>%
        head(200)
      
      # Convertir a EntrezID y realizar ORA
      entrez_ids_mRNAs <- gene_id_to_entrez(top_correlated_mRNAs$gene_id)
      
      if (nrow(entrez_ids_mRNAs) > 0) {
        go_enrich_lncRNA <- clusterProfiler::enrichGO(gene = entrez_ids_mRNAs$ENTREZID,
                                                      OrgDb = org.Hs.eg.db,
                                                      ont = "BP",
                                                      pAdjustMethod = "BH",
                                                      qvalueCutoff = 0.05,
                                                      readable = TRUE)
        
        coexp_output_path <- file.path(results_dir, paste0("09_lncRNA_coexpression_ORA_", lncRNA_id, ".RData"))
        save(go_enrich_lncRNA, file = coexp_output_path)
        message(paste("    - ORA de co-expresión completado. Resultados guardados en:", coexp_output_path))
        
        # --- Generar y guardar plots ORA de co-expresión ---
        if (!is.null(go_enrich_lncRNA) && nrow(go_enrich_lncRNA) > 0) {
          # Dotplot
          dotplot_go <- dotplot(go_enrich_lncRNA, showCategory=15) + ggtitle(paste("Co-expresion ORA:", lncRNA_id))
          ggsave(file.path(results_dir, paste0("09_lncRNA_coexpression_dotplot_", lncRNA_id, ".png")), dotplot_go, width = 10, height = 8)
          
          # Barplot
          barplot_go <- barplot(go_enrich_lncRNA, showCategory=15) + ggtitle(paste("Co-expresion ORA:", lncRNA_id))
          ggsave(file.path(results_dir, paste0("09_lncRNA_coexpression_barplot_", lncRNA_id, ".png")), barplot_go, width = 10, height = 8)
          
          message(paste("    - Gráficos de ORA de co-expresión guardados en:", results_dir))
        }
        
        # --- Guardar resultados ORA de co-expresión en CSV ---
        if (!is.null(go_enrich_lncRNA) && nrow(go_enrich_lncRNA) > 0) {
          results_df <- as.data.frame(go_enrich_lncRNA)
          write.csv(results_df, file = file.path(results_dir, paste0("09_lncRNA_coexpression_ORA_results_", lncRNA_id, ".csv")), row.names = FALSE)
          message(paste("    - Resultados ORA de co-expresión guardados en:", file.path(results_dir, paste0("09_lncRNA_coexpression_ORA_results_", lncRNA_id, ".csv"))))
        }
        
      } else {
        message("    - No se pudieron mapear los mRNAs co-expresados a EntrezID. Saltando ORA.")
      }
    }
  } else {
    message("No se encontraron lncRNAs significativos para analizar.")
  }
} else {
  message("No se encontraron lncRNAs significativos en ninguna de las comparaciones. Saltando el análisis de co-expresión.")
}

message("\nAnálisis de co-expresión completado.")