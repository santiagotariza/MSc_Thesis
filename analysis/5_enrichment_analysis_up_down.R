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
library(enrichplot) # Para DotPlots

message("Librerías cargadas. Iniciando análisis de enriquecimiento separado para genes up y down-regulated...")

# --- FUNCIÓN PARA ENRIQUECIMIENTO ORA (UP/DOWN) POR COMPARACIÓN ---
run_up_down_enrichment <- function(comparison_prefix) {
  message(paste0("\nProcesando enriquecimiento ORA (UP/DOWN) para la comparación: ", comparison_prefix))
  
  # 1. Cargar los DEGs de la comparación específica
  # Asegúrate de que el archivo DEGs_[prefix].csv exista y tenga los IDs de Ensembl como rownames
  degs_current <- read.csv(paste0("DEGs_", comparison_prefix, ".csv"), row.names = 1)
  
  # Extraer los IDs de Ensembl (los rownames del dataframe) y limpiarlos
  original_ensembl_ids <- rownames(degs_current)
  ensembl_ids_clean <- sub("\\..*", "", original_ensembl_ids)
  
  # Añadir los IDs de Ensembl limpios y Entrez IDs al dataframe de DEGs para fácil acceso
  degs_current$ensembl_id_clean <- ensembl_ids_clean
  
  # Mapear a Entrez IDs
  entrez_map <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = degs_current$ensembl_id_clean,
                                      keytype = "ENSEMBL",
                                      columns = "ENTREZID")
  
  # Unir el mapeo con el dataframe de DEGs
  degs_current <- left_join(degs_current, entrez_map, by = c("ensembl_id_clean" = "ENSEMBL")) %>%
    filter(!is.na(ENTREZID)) # Eliminar genes sin Entrez ID
  
  # Separar genes up- y down-regulated
  up_regulated_genes <- degs_current %>% filter(log2FoldChange > 0)
  down_regulated_genes <- degs_current %>% filter(log2FoldChange < 0)
  
  up_entrez_ids <- as.character(up_regulated_genes$ENTREZID)
  down_entrez_ids <- as.character(down_regulated_genes$ENTREZID)
  
  message(paste0("  > Genes UP-REGULATED mapeados a Entrez IDs: ", length(up_entrez_ids)))
  message(paste0("  > Genes DOWN-REGULATED mapeados a Entrez IDs: ", length(down_entrez_ids)))
  
  # --- Enriquecimiento para genes UP-REGULATED ---
  message("  Realizando enriquecimiento para genes UP-REGULATED...")
  if (length(up_entrez_ids) > 0) {
    # KEGG Enrichment para UP-REGULATED
    message("  > KEGG para UP-REGULATED...")
    kegg_up <- enrichKEGG(gene          = up_entrez_ids,
                          organism      = 'hsa',
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 0.05)
    if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
      write.csv(as.data.frame(setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")), paste0("KEGG_enrichment_", comparison_prefix, "_UP.csv"), row.names = FALSE)
      message(paste0("  > Resultados KEGG UP-REGULATED guardados en KEGG_enrichment_", comparison_prefix, "_UP.csv"))
      png(paste0("DotPlot_KEGG_", comparison_prefix, "_UP.png"), width = 900, height = 700, res = 120)
      print(dotplot(kegg_up, showCategory = 15, title = paste0("KEGG Pathway Enrichment (UP-REGULATED, ", comparison_prefix, ")")))
      dev.off()
    } else {
      message(paste0("  > No se encontraron vías KEGG enriquecidas para genes UP-REGULATED en ", comparison_prefix, "."))
    }
    
    # Reactome Enrichment para UP-REGULATED
    message("  > Reactome para UP-REGULATED...")
    reactome_up <- enrichPathway(gene         = up_entrez_ids,
                                 organism     = "human",
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.05,
                                 readable     = TRUE)
    if (!is.null(reactome_up) && nrow(as.data.frame(reactome_up)) > 0) {
      write.csv(as.data.frame(reactome_up), paste0("Reactome_enrichment_", comparison_prefix, "_UP.csv"), row.names = FALSE)
      message(paste0("  > Resultados Reactome UP-REGULATED guardados en Reactome_enrichment_", comparison_prefix, "_UP.csv"))
      png(paste0("DotPlot_Reactome_", comparison_prefix, "_UP.png"), width = 900, height = 700, res = 120)
      print(dotplot(reactome_up, showCategory = 15, title = paste0("Reactome Pathway Enrichment (UP-REGULATED, ", comparison_prefix, ")")))
      dev.off()
    } else {
      message(paste0("  > No se encontraron vías Reactome enriquecidas para genes UP-REGULATED en ", comparison_prefix, "."))
    }
  } else {
    message("  > No hay genes UP-REGULATED para analizar.")
  }
  
  # --- Enriquecimiento para genes DOWN-REGULATED ---
  message("  Realizando enriquecimiento para genes DOWN-REGULATED...")
  if (length(down_entrez_ids) > 0) {
    # KEGG Enrichment para DOWN-REGULATED
    message("  > KEGG para DOWN-REGULATED...")
    kegg_down <- enrichKEGG(gene          = down_entrez_ids,
                            organism      = 'hsa',
                            pAdjustMethod = "BH",
                            qvalueCutoff  = 0.05)
    if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
      write.csv(as.data.frame(setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")), paste0("KEGG_enrichment_", comparison_prefix, "_DOWN.csv"), row.names = FALSE)
      message(paste0("  > Resultados KEGG DOWN-REGULATED guardados en KEGG_enrichment_", comparison_prefix, "_DOWN.csv"))
      png(paste0("DotPlot_KEGG_", comparison_prefix, "_DOWN.png"), width = 900, height = 700, res = 120)
      print(dotplot(kegg_down, showCategory = 15, title = paste0("KEGG Pathway Enrichment (DOWN-REGULATED, ", comparison_prefix, ")")))
      dev.off()
    } else {
      message(paste0("  > No se encontraron vías KEGG enriquecidas para genes DOWN-REGULATED en ", comparison_prefix, "."))
    }
    
    # Reactome Enrichment para DOWN-REGULATED
    message("  > Reactome para DOWN-REGULATED...")
    reactome_down <- enrichPathway(gene          = down_entrez_ids,
                                   organism      = "human",
                                   pAdjustMethod = "BH",
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
    if (!is.null(reactome_down) && nrow(as.data.frame(reactome_down)) > 0) {
      write.csv(as.data.frame(reactome_down), paste0("Reactome_enrichment_", comparison_prefix, "_DOWN.csv"), row.names = FALSE)
      message(paste0("  > Resultados Reactome DOWN-REGULATED guardados en Reactome_enrichment_", comparison_prefix, "_DOWN.csv"))
      png(paste0("DotPlot_Reactome_", comparison_prefix, "_DOWN.png"), width = 900, height = 700, res = 120)
      print(dotplot(reactome_down, showCategory = 15, title = paste0("Reactome Pathway Enrichment (DOWN-REGULATED, ", comparison_prefix, ")")))
      dev.off()
    } else {
      message(paste0("  > No se encontraron vías Reactome enriquecidas para genes DOWN-REGULATED en ", comparison_prefix, "."))
    }
  } else {
    message("  > No hay genes DOWN-REGULATED para analizar.")
  }
  
  message(paste0("Enriquecimiento ORA (UP/DOWN) completado para la comparación: ", comparison_prefix))
}

# --- EJECUTAR PARA CADA COMPARACIÓN ---
# Asegúrate de que los archivos DEGs_WTC_vs_NTC.csv, DEGs_H1KOC_vs_NTC.csv, DEGs_H2KOC_vs_NTC.csv existan
run_up_down_enrichment("WTC_vs_NTC")
run_up_down_enrichment("H1KOC_vs_NTC")
run_up_down_enrichment("H2KOC_vs_NTC")

message("Script 5_enrichment_analysis_up_down.R completado.")