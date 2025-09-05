library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db) # Base de datos de anotación de genes humanos
library(DOSE)
library(GSEABase)
library(ggplot2)
library(stringr)
library(enrichplot) # Paquete para las visualizaciones de enriquecimiento
library(pathview) # Para los mapas de vías KEGG
library(ReactomePA) # Para los análisis de vías de Reactome

# --- CONFIGURACIÓN DE PARÁMETROS ---
# Define los umbrales para filtrar los resultados de la expresión diferencial.
log2FC_threshold <- 1
padj_threshold <- 0.05

# Define el número de lncRNAs principales a analizar para la co-expresión.
top_lncRNAs_to_analyze <- 5

# --- CONFIGURACIÓN DE GRÁFICOS ---
# Activa o desactiva cada tipo de gráfico (TRUE/FALSE)
# Gráficos ORA (GO, KEGG, Reactome)
plot_ora_dotplot <- TRUE
plot_ora_barplot <- TRUE
plot_ora_cnetplot <- TRUE # Red de genes-conceptos
plot_ora_emapplot <- FALSE # Mapa de enriquecimiento (útil para ver solapamiento)
plot_ora_kegg_maps <- TRUE
plot_ora_reactome_plots <- TRUE # Generará dotplot y cnetplot para Reactome

# Gráficos GSEA (GO, KEGG, Reactome)
plot_gsea_gseaplot2 <- TRUE
plot_gsea_dotplot <- TRUE # Muestra los top-pathways de GSEA
plot_gsea_cnetplot <- TRUE
plot_gsea_kegg_maps <- TRUE
plot_gsea_reactome_plots <- TRUE # Generará dotplot y cnetplot para Reactome

# --- FIN DE LA CONFIGURACIÓN ---


# --- Definir directorios de manera robusta ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
}

results_dir <- file.path(script_dir, "../results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

print(results_dir)

# Cargar datos de expresión corregidos y biotipos
load(file.path(results_dir, "01_filtered_data.RData"))
load(file.path(results_dir, "03_dds_vsd_objects.RData"))

# Definir las comparaciones de interés (debe coincidir con el script 06)
comparisons <- c("G3_HMGB1_KO_vs_G1_NTC", "G4_HMGB2_KO_vs_G1_NTC", "G2_WTC_vs_G1_NTC", "G3_HMGB1_KO_vs_G4_HMGB2_KO")

# --- FUNCIONES AUXILIARES ---
# Función para convertir gene_id a EntrezID
gene_id_to_entrez <- function(gene_id_list) {
  bitr(gene_id_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
}

# --- 1. ANÁLISIS DE ENRIQUECIMIENTO (ORA) para mRNA ---
message("\n--- Realizando ORA para genes de mRNA ---")

for (comp_name in comparisons) {
  message(paste("  - Procesando ORA para:", comp_name))
  
  de_results_path <- file.path(results_dir, paste0("06_mRNA_", comp_name, "_DE_results.csv"))
  if (!file.exists(de_results_path)) {
    warning(paste("Archivo no encontrado:", de_results_path, ". Saltando ORA para esta comparación."))
    next
  }
  
  de_results <- read.csv(de_results_path, stringsAsFactors = FALSE)
  significant_genes <- de_results$gene_id
  entrez_ids <- gene_id_to_entrez(significant_genes)
  
  if (nrow(entrez_ids) > 0) {
    go_enrich <- clusterProfiler::enrichGO(gene = entrez_ids$ENTREZID,
                                           OrgDb = org.Hs.eg.db,
                                           ont = "BP",
                                           pAdjustMethod = "BH",
                                           qvalueCutoff = 0.05,
                                           readable = TRUE)
    
    kegg_enrich <- clusterProfiler::enrichKEGG(gene = entrez_ids$ENTREZID,
                                               organism = 'hsa',
                                               pAdjustMethod = "BH",
                                               qvalueCutoff = 0.05)
    
    # Análisis ORA para Reactome
    reactome_enrich <- ReactomePA::enrichPathway(gene = entrez_ids$ENTREZID,
                                                 pAdjustMethod = "BH",
                                                 qvalueCutoff = 0.05,
                                                 readable = TRUE)
    
    ora_output_path <- file.path(results_dir, paste0("07_mRNA_ORA_", comp_name, ".RData"))
    save(go_enrich, kegg_enrich, reactome_enrich, file = ora_output_path)
    message(paste("    - ORA completado. Resultados guardados en:", ora_output_path))
    
    # --- Generar y guardar plots ORA basados en la configuración ---
    if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
      if (plot_ora_dotplot) {
        dotplot_go <- enrichplot::dotplot(go_enrich, showCategory=15) + ggtitle(paste("ORA Dotplot:", comp_name))
        ggsave(file.path(results_dir, paste0("07_mRNA_GO_dotplot_", comp_name, ".png")), dotplot_go, width = 10, height = 8)
      }
      if (plot_ora_barplot) {
        barplot_go <- clusterProfiler::barplot(go_enrich, showCategory=15) + ggtitle(paste("ORA Barplot:", comp_name))
        ggsave(file.path(results_dir, paste0("07_mRNA_GO_barplot_", comp_name, ".png")), barplot_go, width = 10, height = 8)
      }
      if (plot_ora_cnetplot) {
        cnetplot_go <- enrichplot::cnetplot(go_enrich, foldChange = de_results$log2FoldChange, node_label="category") + ggtitle(paste("ORA Cnetplot:", comp_name))
        ggsave(file.path(results_dir, paste0("07_mRNA_GO_cnetplot_", comp_name, ".png")), cnetplot_go, width = 10, height = 8)
      }
      if (plot_ora_emapplot) {
        if (nrow(go_enrich) > 1) {
          emapplot_go <- enrichplot::emapplot(go_enrich) + ggtitle(paste("ORA Emapplot:", comp_name))
          ggsave(file.path(results_dir, paste0("07_mRNA_GO_emapplot_", comp_name, ".png")), emapplot_go, width = 10, height = 8)
        } else {
          message("    - No hay suficientes vías para generar un Emapplot. Saltando.")
        }
      }
    }
    
    if (plot_ora_reactome_plots && !is.null(reactome_enrich) && nrow(reactome_enrich) > 0) {
      if (plot_ora_dotplot) {
        dotplot_reactome <- enrichplot::dotplot(reactome_enrich, showCategory=15) + ggtitle(paste("ORA Reactome Dotplot:", comp_name))
        ggsave(file.path(results_dir, paste0("07_mRNA_Reactome_dotplot_", comp_name, ".png")), dotplot_reactome, width = 10, height = 8)
      }
      if (plot_ora_cnetplot) {
        cnetplot_reactome <- enrichplot::cnetplot(reactome_enrich, foldChange = de_results$log2FoldChange, node_label="category") + ggtitle(paste("ORA Reactome Cnetplot:", comp_name))
        ggsave(file.path(results_dir, paste0("07_mRNA_Reactome_cnetplot_", comp_name, ".png")), cnetplot_reactome, width = 10, height = 8)
      }
    }
    
    # --- Generar mapas de vías KEGG ---
    if (plot_ora_kegg_maps && !is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
      logFC_vector_ora <- de_results$log2FoldChange
      names(logFC_vector_ora) <- de_results$gene_id
      
      top_kegg_pathways <- head(kegg_enrich$ID, 3)
      for (pathway_id in top_kegg_pathways) {
        tryCatch({
          pathview::pathview(gene.data = logFC_vector_ora,
                             pathway.id = pathway_id,
                             species = "hsa",
                             gene.idtype = "ENSEMBL",
                             kegg.native = TRUE,
                             out.suffix = paste0("ORA_mRNA_KEGG_", comp_name),
                             limit = list(gene=2, cpd=1),
                             out.dir = results_dir) # **CORRECCIÓN AQUÍ**
          message(paste("    - Mapa de vía KEGG generado para:", pathway_id))
        }, error = function(e) {
          warning(paste("Error al generar el mapa de vía para", pathway_id, ":", e$message))
        })
      }
    }
  } else {
    message("    - No se pudieron mapear genes significativos a EntrezID. Saltando ORA.")
  }
}

# --- 2. ANÁLISIS DE ENRIQUECIMIENTO (GSEA) para mRNA ---
message("\n--- Realizando GSEA para genes de mRNA ---")

for (comp_name in comparisons) {
  message(paste("  - Procesando GSEA para:", comp_name))
  
  de_results_path <- file.path(results_dir, paste0("06_mRNA_", comp_name, "_DE_results.csv"))
  if (!file.exists(de_results_path)) {
    warning(paste("Archivo no encontrado:", de_results_path, ". Saltando GSEA para esta comparación."))
    next
  }
  
  de_results <- read.csv(de_results_path, stringsAsFactors = FALSE)
  dds_path <- file.path(results_dir, "03_dds_vsd_objects.RData")
  if (!file.exists(dds_path)) {
    warning("Objeto DESeqDataSet no encontrado. Saltando GSEA.")
    next
  }
  load(dds_path)
  dds <- mRNA_obj$dds
  results_dds <- results(dds, contrast = c("condition", str_split(comp_name, "_vs_")[[1]]), alpha = padj_threshold)
  
  gene_list <- results_dds$stat
  names(gene_list) <- rownames(results_dds)
  gene_list <- gene_list[!is.na(gene_list)]
  
  gene_list_df <- data.frame(ENSEMBL = names(gene_list), stat = gene_list)
  gene_list_df <- gene_list_df %>% 
    left_join(bitr(unique(gene_list_df$ENSEMBL), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db"), by = "ENSEMBL") %>%
    filter(!is.na(ENTREZID))
  
  gene_list_filtered <- gene_list_df %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(stat), n = 1) %>%
    ungroup() %>%
    dplyr::select(ENTREZID, stat)
  
  final_gene_list <- gene_list_filtered$stat
  names(final_gene_list) <- gene_list_filtered$ENTREZID
  final_gene_list <- sort(final_gene_list, decreasing = TRUE)
  
  if (length(final_gene_list) > 0) {
    gsea_go <- clusterProfiler::gseGO(geneList = final_gene_list, 
                                      OrgDb = org.Hs.eg.db,
                                      ont = "BP",
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      pvalueCutoff = 0.05,
                                      verbose = FALSE)
    
    gsea_kegg <- clusterProfiler::gseKEGG(geneList = final_gene_list,
                                          organism = 'hsa',
                                          minGSSize = 10,
                                          maxGSSize = 500,
                                          pvalueCutoff = 0.05,
                                          verbose = FALSE)
    
    # Análisis GSEA para Reactome
    gsea_reactome <- ReactomePA::gsePathway(geneList = final_gene_list,
                                            pvalueCutoff = 0.05,
                                            verbose = FALSE)
    
    gsea_output_path <- file.path(results_dir, paste0("07_mRNA_GSEA_", comp_name, ".RData"))
    save(gsea_go, gsea_kegg, gsea_reactome, file = gsea_output_path)
    message(paste("    - GSEA completado. Resultados guardados en:", gsea_output_path))
    
    # --- Generar y guardar plots GSEA basados en la configuración ---
    if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
      if (plot_gsea_gseaplot2) {
        gsea_plot <- enrichplot::gseaplot2(gsea_go, geneSetID = 1:min(5, nrow(gsea_go)), title = paste("GSEA - Top Pathways for mRNA:", comp_name))
        ggsave(file.path(results_dir, paste0("07_mRNA_GSEA_gseaplot2_", comp_name, ".png")), gsea_plot, width = 10, height = 8)
      }
      if (plot_gsea_dotplot) {
        dotplot_gsea <- enrichplot::dotplot(gsea_go, showCategory = 15) + ggtitle(paste("GSEA Dotplot:", comp_name))
        ggsave(file.path(results_dir, paste0("07_mRNA_GSEA_dotplot_", comp_name, ".png")), dotplot_gsea, width = 10, height = 8)
      }
      if (plot_gsea_cnetplot) {
        cnetplot_gsea <- enrichplot::cnetplot(gsea_go, foldChange = final_gene_list, node_label="category") + ggtitle(paste("GSEA Cnetplot:", comp_name))
        ggsave(file.path(results_dir, paste0("07_mRNA_GSEA_cnetplot_", comp_name, ".png")), cnetplot_gsea, width = 10, height = 8)
      }
    }
    
    if (plot_gsea_reactome_plots && !is.null(gsea_reactome) && nrow(gsea_reactome) > 0) {
      if (plot_gsea_dotplot) {
        dotplot_reactome_gsea <- enrichplot::dotplot(gsea_reactome, showCategory = 15) + ggtitle(paste("GSEA Reactome Dotplot:", comp_name))
        ggsave(file.path(results_dir, paste0("07_mRNA_GSEA_Reactome_dotplot_", comp_name, ".png")), dotplot_reactome_gsea, width = 10, height = 8)
      }
      if (plot_gsea_cnetplot) {
        cnetplot_reactome_gsea <- enrichplot::cnetplot(gsea_reactome, foldChange = final_gene_list, node_label="category") + ggtitle(paste("GSEA Reactome Cnetplot:", comp_name))
        ggsave(file.path(results_dir, paste0("07_mRNA_GSEA_Reactome_cnetplot_", comp_name, ".png")), cnetplot_reactome_gsea, width = 10, height = 8)
      }
    }
    
    # --- Generar mapas de vías KEGG ---
    if (plot_gsea_kegg_maps && !is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
      top_kegg_pathways <- head(gsea_kegg$ID, 3)
      for (pathway_id in top_kegg_pathways) {
        tryCatch({
          pathview::pathview(gene.data = final_gene_list,
                             pathway.id = pathway_id,
                             species = "hsa",
                             gene.idtype = "ENTREZID",
                             kegg.native = TRUE,
                             out.suffix = paste0("GSEA_mRNA_KEGG_", comp_name),
                             out.dir = results_dir) # **CORRECCIÓN AQUÍ**
          message(paste("    - Mapa de vía KEGG generado para:", pathway_id))
        }, error = function(e) {
          warning(paste("Error al generar el mapa de vía para", pathway_id, ":", e$message))
        })
      }
    }
  } else {
    message("    - No se pudieron mapear genes para GSEA. Saltando análisis.")
  }
}


# --- 3. ANÁLISIS DE CO-EXPRESIÓN para lncRNA ---
message("\n--- Realizando análisis de co-expresión para lncRNAs ---")

# Cargar los datos de lncRNA y mRNA vsd
vsd_mRNA <- mRNA_obj$vsd
vsd_lncRNA <- lncRNA_obj$vsd

# Obtener una lista maestra de lncRNAs significativos de todas las comparaciones
all_lncRNAs_significant <- data.frame()
lncRNA_comparisons <- c("G3_HMGB1_KO_vs_G1_NTC", "G4_HMGB2_KO_vs_G1_NTC", "G2_WTC_vs_G1_NTC", "G3_HMGB1_KO_vs_G4_HMGB2_KO")

for (comp_name in lncRNA_comparisons) {
  lncRNA_de_results_path <- file.path(results_dir, paste0("06_lncRNA_", comp_name, "_DE_results.csv"))
  
  if (file.exists(lncRNA_de_results_path)) {
    de_results <- read.csv(lncRNA_de_results_path, stringsAsFactors = FALSE)
    if (nrow(de_results) > 0) {
      de_results$comparison <- comp_name # Añadir la columna de comparación
      all_lncRNAs_significant <- bind_rows(all_lncRNAs_significant, de_results)
    }
  }
}

if (nrow(all_lncRNAs_significant) > 0) {
  # Eliminar duplicados, manteniendo el resultado con el padj más bajo
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
        head(200) # Se toman los 200 más correlacionados
      
      # Convertir a EntrezID y realizar ORA
      entrez_ids_mRNAs <- gene_id_to_entrez(top_correlated_mRNAs$gene_id)
      
      if (nrow(entrez_ids_mRNAs) > 0) {
        # ORA para los genes co-expresados
        go_enrich_lncRNA <- clusterProfiler::enrichGO(gene = entrez_ids_mRNAs$ENTREZID,
                                                      OrgDb = org.Hs.eg.db,
                                                      ont = "BP",
                                                      pAdjustMethod = "BH",
                                                      qvalueCutoff = 0.05,
                                                      readable = TRUE)
        
        # Guardar resultados del ORA de co-expresión
        coexp_output_path <- file.path(results_dir, paste0("07_lncRNA_coexpression_ORA_", lncRNA_id, ".RData"))
        save(go_enrich_lncRNA, file = coexp_output_path)
        message(paste("    - ORA de co-expresión completado. Resultados guardados en:", coexp_output_path))
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

message("\nAnálisis funcional y de co-expresión completado. Revisa la carpeta 'results' para los resultados y los gráficos.")