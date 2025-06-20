# Instalar librerías si no están instaladas
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("fgsea", quietly = TRUE))
  BiocManager::install("fgsea")
if (!requireNamespace("msigdbr", quietly = TRUE))
  install.packages("msigdbr") # Para descargar conjuntos de genes
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

library(fgsea)
library(msigdbr)
library(dplyr) # Cargar dplyr explícitamente
library(ggplot2)
library(org.Hs.eg.db) # Cargar para mapeo de IDs

message("Librerías para GSEA cargadas. Iniciando análisis GSEA...")

# --- Función para preparar los datos y ejecutar GSEA ---
run_gsea <- function(file_path, comparison_name) {
  message(paste0("\nProcesando archivo para GSEA: ", file_path))
  
  # 1. Cargar los resultados completos de DESeq2
  df_res <- read.csv(file_path, row.names = 1)
  
  # Preparar la lista de clasificación (rank list)
  # Necesitamos un vector numérico con LFC, nombrado por Entrez ID.
  # Los rownames son los Ensembl ID (con versiones, ej. ENSG0...
  
  # Limpiar Ensembl IDs de la versión
  original_ensembl_ids <- rownames(df_res)
  ensembl_ids_clean <- sub("\\..*", "", original_ensembl_ids) # Eliminar ".1", ".2", etc.
  
  # Mapear Ensembl IDs a Entrez IDs
  gene_list <- df_res$log2FoldChange
  names(gene_list) <- ensembl_ids_clean
  
  # Eliminar duplicados si los hay (tomar el primero o el promedio, aquí tomamos el promedio)
  # (Se usa un enfoque diferente aquí que en los scripts ORA para manejar LFC de duplicados)
  gene_list_df <- data.frame(ENSEMBL = names(gene_list), log2FoldChange = gene_list)
  
  # Mapear a Entrez IDs
  entrez_map <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = names(gene_list),
                                      keytype = "ENSEMBL",
                                      columns = "ENTREZID")
  
  gene_list_df <- left_join(gene_list_df, entrez_map, by = c("ENSEMBL" = "ENSEMBL")) %>%
    filter(!is.na(ENTREZID)) %>% # Eliminar genes sin Entrez ID
    group_by(ENTREZID) %>%
    summarise(log2FoldChange = mean(log2FoldChange)) # Promediar si un Entrez ID mapea a múltiples Ensembl
  
  # Convertir a un vector nombrado para fgsea
  ranked_genes <- setNames(gene_list_df$log2FoldChange, gene_list_df$ENTREZID)
  
  # Ordenar la lista de genes por LFC (descendente)
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  message(paste0("  Número de genes rankeados para GSEA: ", length(ranked_genes)))
  
  # 2. Cargar conjuntos de genes (Pathways) de MSigDB para KEGG
  m_t2g_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY") %>%
    dplyr::select(gs_name, entrez_gene)
  
  # Convertir a una lista de pathways
  pathways_kegg <- split(m_t2g_kegg$entrez_gene, f = m_t2g_kegg$gs_name)
  
  # Realizar GSEA para KEGG
  message("  Realizando GSEA para vías KEGG...")
  gsea_kegg_res <- fgsea(pathways = pathways_kegg,
                         stats = ranked_genes,
                         minSize = 10,
                         maxSize = 500,
                         nPermSimple = 10000) # Aumentado para más robustez
  
  # Filtrar y guardar resultados de GSEA para KEGG
  gsea_kegg_res_filtered <- as.data.frame(gsea_kegg_res) %>%
    filter(padj < 0.05) %>%
    dplyr::select(pathway, NES, pval, padj, ES, size) %>% # Excluir 'leadingEdge' y 'pathways' (list columns)
    arrange(NES)
  
  write.csv(gsea_kegg_res_filtered, paste0("GSEA_KEGG_", comparison_name, ".csv"), row.names = FALSE)
  message(paste0("  > Resultados GSEA KEGG guardados en GSEA_KEGG_", comparison_name, ".csv"))
  
  if (nrow(gsea_kegg_res_filtered) > 0) {
    message("  Generando DotPlot para GSEA KEGG...")
    png(paste0("DotPlot_GSEA_KEGG_", comparison_name, ".png"), width = 1200, height = 700, res = 120) # Aumenta el ancho
    print(ggplot(gsea_kegg_res_filtered %>% dplyr::top_n(15, wt = abs(NES)), aes(x = NES, y = reorder(pathway, NES), color = padj, size = size)) +
            geom_point() +
            scale_color_viridis_c(trans = "log10", direction = -1) +
            labs(title = paste0("GSEA KEGG Pathways (", comparison_name, ")"),
                 x = "Normalized Enrichment Score (NES)", y = "Pathway") +
            theme_minimal() +
            theme(legend.position = "right",
                  plot.title = element_text(hjust = 0.5),
                  axis.text.y = element_text(size = 7)))
    dev.off()
    message(paste0("  > DotPlot GSEA KEGG guardado en DotPlot_GSEA_KEGG_", comparison_name, ".png"))
  } else {
    message("  > No se encontraron vías KEGG significativas en GSEA para ", comparison_name, ".")
  }
  
  # 3. Cargar conjuntos de genes para Reactome
  m_t2g_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
    dplyr::select(gs_name, entrez_gene)
  
  pathways_reactome <- split(m_t2g_reactome$entrez_gene, f = m_t2g_reactome$gs_name)
  
  # Realizar GSEA para Reactome
  message("  Realizando GSEA para vías Reactome...")
  gsea_reactome_res <- fgsea(pathways = pathways_reactome,
                             stats = ranked_genes,
                             minSize = 10,
                             maxSize = 500,
                             nPermSimple = 10000)
  
  # Filtrar y guardar resultados de GSEA para Reactome
  gsea_reactome_res_filtered <- as.data.frame(gsea_reactome_res) %>%
    filter(padj < 0.05) %>%
    dplyr::select(pathway, NES, pval, padj, ES, size) %>% # Excluir 'leadingEdge' y 'pathways' (list columns)
    arrange(NES)
  
  write.csv(gsea_reactome_res_filtered, paste0("GSEA_Reactome_", comparison_name, ".csv"), row.names = FALSE)
  message(paste0("  > Resultados GSEA Reactome guardados en GSEA_Reactome_", comparison_name, ".csv"))
  
  if (nrow(gsea_reactome_res_filtered) > 0) {
    message("  Generando DotPlot para GSEA Reactome...")
    png(paste0("DotPlot_GSEA_Reactome_", comparison_name, ".png"), width = 1200, height = 700, res = 120) # Aumenta el ancho
    print(ggplot(gsea_reactome_res_filtered %>% dplyr::top_n(15, wt = abs(NES)), aes(x = NES, y = reorder(pathway, NES), color = padj, size = size)) +
            geom_point() +
            scale_color_viridis_c(trans = "log10", direction = -1) +
            labs(title = paste0("GSEA Reactome Pathways (", comparison_name, ")"),
                 x = "Normalized Enrichment Score (NES)", y = "Pathway") +
            theme_minimal() +
            theme(legend.position = "right",
                  plot.title = element_text(hjust = 0.5),
                  axis.text.y = element_text(size = 7)))
    dev.off()
    message(paste0("  > DotPlot GSEA Reactome guardado en DotPlot_GSEA_Reactome_", comparison_name, ".png"))
  } else {
    message("  > No se encontraron vías Reactome significativas en GSEA para ", comparison_name, ".")
  }
}

# --- Ejecutar GSEA para cada comparación (¡AJUSTADAS LAS COMPARACIONES AQUÍ!) ---
# Asegúrate de que estos archivos existan en tu directorio de trabajo
run_gsea("DESeq2_results_full_WTC_vs_NTC.csv", "WTC_vs_NTC")
run_gsea("DESeq2_results_full_H1KOC_vs_NTC.csv", "H1KOC_vs_NTC")
run_gsea("DESeq2_results_full_H2KOC_vs_NTC.csv", "H2KOC_vs_NTC")

message("Script 7_GSEA_analysis.R completado.")