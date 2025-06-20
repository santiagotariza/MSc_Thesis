# Instalar librerías si no están instaladas
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA")
if (!requireNamespace("DOSE", quietly = TRUE))
  BiocManager::install("DOSE") # Dependencia de clusterProfiler

# Cargar librerías
library(clusterProfiler)
library(org.Hs.eg.db) # Base de datos de anotación para Homo sapiens
library(ReactomePA)
library(dplyr) # Para manipulación de datos

message("Librerías cargadas. Iniciando análisis de enriquecimiento...")

# --- FUNCIÓN PARA ENRIQUECIMIENTO ORA POR COMPARACIÓN ---
run_ora_enrichment <- function(comparison_prefix) {
  message(paste0("\nProcesando enriquecimiento ORA para la comparación: ", comparison_prefix))
  
  # 1. Cargar los DEGs significativos de la comparación específica
  # Asegúrate de que el archivo DEGs_[prefix].csv exista
  degs_current <- read.csv(paste0("DEGs_", comparison_prefix, ".csv"), row.names = 1)
  
  # Extraer los IDs de Ensembl (los rownames del dataframe) y limpiarlos
  original_ensembl_ids <- rownames(degs_current)
  ensembl_ids_clean <- sub("\\..*", "", original_ensembl_ids) # Eliminar ".1", ".2", etc.
  
  # Mapear IDs de Ensembl a Entrez IDs
  # Usar select para mapeo 1:1, si hay múltiples Entrez por Ensembl, se seleccionará el primero
  entrez_map <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = ensembl_ids_clean,
                                      keytype = "ENSEMBL",
                                      columns = "ENTREZID")
  
  # Filtrar los genes que no se pudieron mapear a Entrez ID
  entrez_ids <- na.omit(entrez_map$ENTREZID)
  
  message(paste0("  Número de DEGs mapeados a Entrez IDs para ", comparison_prefix, ": ", length(entrez_ids)))
  
  if (length(entrez_ids) == 0) {
    message(paste0("  No hay Entrez IDs válidos para ", comparison_prefix, ". Saltando el análisis de enriquecimiento."))
    return(NULL)
  }
  
  # 2. Análisis de Enriquecimiento de GO (Biological Process)
  message("  Realizando análisis de enriquecimiento GO (BP)...")
  go_results_bp <- enrichGO(gene         = as.character(entrez_ids), # Convertir a character
                            OrgDb        = org.Hs.eg.db,
                            keyType      = 'ENTREZID',
                            ont          = "BP", # Biological Process
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable     = TRUE) # Para obtener Gene Symbols
  
  if (!is.null(go_results_bp) && nrow(as.data.frame(go_results_bp)) > 0) {
    go_bp_df <- as.data.frame(go_results_bp)
    write.csv(go_bp_df, paste0("GO_enrichment_", comparison_prefix, "_BP.csv"), row.names = FALSE)
    message(paste0("  > Resultados GO (BP) guardados en GO_enrichment_", comparison_prefix, "_BP.csv"))
  } else {
    message(paste0("  > No se encontraron vías GO (BP) enriquecidas para ", comparison_prefix, "."))
  }
  
  # 3. Análisis de Enriquecimiento de GO (Molecular Function)
  message("  Realizando análisis de enriquecimiento GO (MF)...")
  go_results_mf <- enrichGO(gene         = as.character(entrez_ids),
                            OrgDb        = org.Hs.eg.db,
                            keyType      = 'ENTREZID',
                            ont          = "MF", # Molecular Function
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable     = TRUE)
  
  if (!is.null(go_results_mf) && nrow(as.data.frame(go_results_mf)) > 0) {
    go_mf_df <- as.data.frame(go_results_mf)
    write.csv(go_mf_df, paste0("GO_enrichment_", comparison_prefix, "_MF.csv"), row.names = FALSE)
    message(paste0("  > Resultados GO (MF) guardados en GO_enrichment_", comparison_prefix, "_MF.csv"))
  } else {
    message(paste0("  > No se encontraron vías GO (MF) enriquecidas para ", comparison_prefix, "."))
  }
  
  # 4. Análisis de Enriquecimiento de GO (Cellular Component)
  message("  Realizando análisis de enriquecimiento GO (CC)...")
  go_results_cc <- enrichGO(gene         = as.character(entrez_ids),
                            OrgDb        = org.Hs.eg.db,
                            keyType      = 'ENTREZID',
                            ont          = "CC", # Cellular Component
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable     = TRUE)
  
  if (!is.null(go_results_cc) && nrow(as.data.frame(go_results_cc)) > 0) {
    go_cc_df <- as.data.frame(go_results_cc)
    write.csv(go_cc_df, paste0("GO_enrichment_", comparison_prefix, "_CC.csv"), row.names = FALSE)
    message(paste0("  > Resultados GO (CC) guardados en GO_enrichment_", comparison_prefix, "_CC.csv"))
  } else {
    message(paste0("  > No se encontraron vías GO (CC) enriquecidas para ", comparison_prefix, "."))
  }
  
  # 5. Análisis de Enriquecimiento de Vías KEGG
  message("  Realizando análisis de enriquecimiento KEGG...")
  kegg_results <- enrichKEGG(gene         = as.character(entrez_ids),
                             organism     = 'hsa', # Homo sapiens
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05)
  
  if (!is.null(kegg_results) && nrow(as.data.frame(kegg_results)) > 0) {
    kegg_df <- as.data.frame(setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
    write.csv(kegg_df, paste0("KEGG_enrichment_", comparison_prefix, ".csv"), row.names = FALSE)
    message(paste0("  > Resultados KEGG guardados en KEGG_enrichment_", comparison_prefix, ".csv"))
  } else {
    message(paste0("  > No se encontraron vías KEGG enriquecidas para ", comparison_prefix, "."))
  }
  
  # 6. Análisis de Enriquecimiento de Vías Reactome
  message("  Realizando análisis de enriquecimiento Reactome...")
  reactome_results <- enrichPathway(gene         = as.character(entrez_ids),
                                    organism     = "human",
                                    pAdjustMethod = "BH",
                                    qvalueCutoff = 0.05,
                                    readable     = TRUE) # Para Gene Symbols
  
  if (!is.null(reactome_results) && nrow(as.data.frame(reactome_results)) > 0) {
    reactome_df <- as.data.frame(reactome_results)
    write.csv(reactome_df, paste0("Reactome_enrichment_", comparison_prefix, ".csv"), row.names = FALSE)
    message(paste0("  > Resultados Reactome guardados en Reactome_enrichment_", comparison_prefix, ".csv"))
  } else {
    message(paste0("  > No se encontraron vías Reactome enriquecidas para ", comparison_prefix, "."))
  }
  
  message(paste0("Enriquecimiento ORA completado para la comparación: ", comparison_prefix))
}

# --- EJECUTAR PARA CADA COMPARACIÓN ---
# Asegúrate de que los archivos DEGs_WTC_vs_NTC.csv, DEGs_H1KOC_vs_NTC.csv, DEGs_H2KOC_vs_NTC.csv existan
run_ora_enrichment("WTC_vs_NTC")
run_ora_enrichment("H1KOC_vs_NTC")
run_ora_enrichment("H2KOC_vs_NTC")

message("Script 3_enrichment_analysis.R completado.")