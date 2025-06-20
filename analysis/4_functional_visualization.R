# Instalar librerías si no están instaladas
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("enrichplot", quietly = TRUE))
  BiocManager::install("enrichplot")
if (!requireNamespace("pathview", quietly = TRUE)) # Para visualizar vías KEGG
  BiocManager::install("pathview")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA")
if (!requireNamespace("dplyr", quietly = TRUE)) # Para manipulación de datos
  install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) # Para gráficos base
  install.packages("ggplot2")
if (!requireNamespace("tibble", quietly = TRUE)) # Para rownames_to_column
  install.packages("tibble")

# Cargar librerías
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot) # Para visualizaciones de enriquecimiento
library(ggplot2) # Para gráficos base
library(dplyr)
library(pathview) # Opcional: para visualizar mapas de vías KEGG detallados
library(tibble) # Para rownames_to_column

message("Iniciando visualización de resultados de enriquecimiento...")

# --- FUNCIÓN PARA GENERAR VISUALIZACIONES POR COMPARACIÓN ---
generate_enrichment_visualizations <- function(comparison_prefix) {
  message(paste0("\nGenerando visualizaciones para la comparación: ", comparison_prefix))
  
  # --- A) Cargar los resultados de expresión diferencial completos para LFC ---
  # Esto es necesario para cnetplot y pathview, que requieren el log2FoldChange de todos los genes.
  full_deg_results_path <- paste0("DESeq2_results_full_", comparison_prefix, ".csv")
  if (file.exists(full_deg_results_path)) {
    full_deg_df <- read.csv(full_deg_results_path)
    
    # Preparar el vector de log2FoldChange con Entrez IDs para cnetplot y pathview
    # Limpiar Ensembl IDs de la versión
    full_deg_df$ensembl_id_clean <- sub("\\..*", "", full_deg_df$GeneID)
    
    # Mapear Ensembl IDs a Entrez IDs
    entrez_map <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = full_deg_df$ensembl_id_clean,
                                        keytype = "ENSEMBL",
                                        columns = "ENTREZID")
    # Unir el mapeo y filtrar NA, luego agrupar por EntrezID y tomar el promedio del LFC
    fc_vector_df <- left_join(full_deg_df, entrez_map, by = c("ensembl_id_clean" = "ENSEMBL")) %>%
      filter(!is.na(ENTREZID)) %>%
      group_by(ENTREZID) %>%
      summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE)) # na.rm=TRUE por si hay NAs en LFC
    
    fc_vector_entrez <- setNames(fc_vector_df$log2FoldChange, fc_vector_df$ENTREZID)
    message(paste0("  > Se preparó el vector de LFC con Entrez IDs para ", comparison_prefix))
    
  } else {
    message(paste0("  ADVERTENCIA: No se encontró el archivo de resultados completos de DEGs: ", full_deg_results_path))
    message("  > No se generarán cnetplots ni pathview para esta comparación.")
    fc_vector_entrez <- NULL
  }
  
  # --- B) Cargar resultados de enriquecimiento ORA (general) ---
  go_bp_path <- paste0("GO_enrichment_", comparison_prefix, "_BP.csv")
  kegg_path <- paste0("KEGG_enrichment_", comparison_prefix, ".csv")
  reactome_path <- paste0("Reactome_enrichment_", comparison_prefix, ".csv")
  
  go_results_bp <- NULL
  if (file.exists(go_bp_path)) {
    go_results_bp <- read.csv(go_bp_path)
    # Convertir a objeto 'enrichResult' si se cargó desde CSV para poder usar enrichplot
    # Asegurarse de que el dataframe no esté vacío antes de convertir
    if(nrow(go_results_bp) > 0) {
      go_results_bp <- clusterProfiler::build_enrichResult(go_results_bp)
    } else {
      go_results_bp <- NULL # Establecer a NULL si está vacío
    }
  }
  
  kegg_results <- NULL
  if (file.exists(kegg_path)) {
    kegg_results <- read.csv(kegg_path)
    if(nrow(kegg_results) > 0) {
      kegg_results <- clusterProfiler::build_enrichResult(kegg_results)
    } else {
      kegg_results <- NULL
    }
  }
  
  reactome_results <- NULL
  if (file.exists(reactome_path)) {
    reactome_results <- read.csv(reactome_path)
    if(nrow(reactome_results) > 0) {
      reactome_results <- clusterProfiler::build_enrichResult(reactome_results)
    } else {
      reactome_results <- NULL
    }
  }
  
  # --- C) Generar Plots para enriquecimiento ORA (general) ---
  
  # DotPlot para GO (BP)
  if (!is.null(go_results_bp) && nrow(as.data.frame(go_results_bp)) > 0) {
    message("  Generando DotPlot para GO (BP)...")
    png(paste0("DotPlot_GO_", comparison_prefix, "_BP.png"), width = 900, height = 700, res = 120)
    print(dotplot(go_results_bp, showCategory = 15, title = paste0("GO Biological Process Enrichment (", comparison_prefix, ")")))
    dev.off()
  } else {
    message(paste0("  > No se encontraron resultados GO (BP) para plotear para ", comparison_prefix, "."))
  }
  
  # DotPlot para KEGG
  if (!is.null(kegg_results) && nrow(as.data.frame(kegg_results)) > 0) {
    message("  Generando DotPlot para KEGG... (General)")
    png(paste0("DotPlot_KEGG_", comparison_prefix, ".png"), width = 900, height = 700, res = 120)
    print(dotplot(kegg_results, showCategory = 15, title = paste0("KEGG Pathway Enrichment (", comparison_prefix, ")")))
    dev.off()
    
    # Cnetplot para KEGG (requiere fc_vector_entrez)
    if (!is.null(fc_vector_entrez)) {
      message("  Generando CnetPlot para KEGG...")
      tryCatch({
        png(paste0("CnetPlot_KEGG_", comparison_prefix, ".png"), width = 1000, height = 900, res = 120)
        # Limitar a las 5 vías más significativas si hay muchas para que el plot no sea ilegible
        kegg_for_cnet <- kegg_results
        if (nrow(as.data.frame(kegg_results)) > 5) {
          kegg_for_cnet_df <- as.data.frame(kegg_results)
          kegg_for_cnet_df <- top_n(kegg_for_cnet_df, 5, wt = -p.adjust)
          kegg_for_cnet <- clusterProfiler::build_enrichResult(kegg_for_cnet_df) # Convertir de nuevo a enrichResult
        }
        
        # Necesitamos que los IDs de los genes en el objeto kegg_for_cnet sean Entrez IDs para cnetplot
        # y que coincidan con los nombres en fc_vector_entrez.
        # clusterProfiler::setReadable transforma a SYMBOL, pero cnetplot usa los IDs originales internamente
        # y mapea los foldChange con los nombres del vector.
        # Asegurarse que el campo 'geneID' en kegg_for_cnet es EntrezID
        
        # Una forma más robusta de asegurar que cnetplot funciona bien:
        # Convertir a df, filtrar top pathways, asegurarse de que GeneID contiene Entrez
        # Reconstruir el objeto enrichResult
        kegg_df_cnet_prep <- as.data.frame(kegg_for_cnet)
        
        # Convertir la columna 'geneID' de string con '/' a una lista de Entrez IDs
        kegg_df_cnet_prep$geneID <- lapply(strsplit(kegg_df_cnet_prep$geneID, "/"), as.character)
        
        # Reconstruir un objeto enrichResult con los EntrezID para cnetplot
        # Esto es un poco hacky porque build_enrichResult espera el formato original
        # Una alternativa es asegurar que el objeto kegg_results *antes* del setReadable
        # tenga los EntrezID y usar ese para cnetplot.
        
        # Para simplificar y dado que el fc_vector ya tiene EntrezID,
        # aseguremos que el objeto de clusterProfiler (kegg_results)
        # tiene los geneIDs como EntrezIDs (que debería ser el caso si no se usó setReadable antes de plotear)
        
        # Si kegg_results ya fue `setReadable` con Gene Symbols, cnetplot podría tener problemas.
        # Volver a cargar o asegurarse de usar el objeto `kegg_results` original (con Entrez IDs).
        
        # Para cnetplot, el objeto de enriquecimiento debe tener los Entrez IDs en la columna 'geneID'
        # Si clusterProfiler::setReadable ya fue aplicado, esa columna podría contener Gene Symbols.
        # Hay que tener cuidado aquí.
        
        # Recomiendo que el objeto que pasas a cnetplot sea el objeto enrichResult original (con Entrez IDs)
        # antes de convertir a readable, si es que lo has hecho.
        
        # Si se usa `setReadable` para convertir a símbolos para los archivos CSV,
        # entonces para cnetplot se necesitaría volver a cargar los datos originales
        # o guardar una versión sin `setReadable` para visualizaciones.
        
        # Asumiendo que `kegg_results` (el objeto clusterProfiler) aún contiene Entrez IDs en su forma interna:
        print(cnetplot(kegg_results,
                       foldChange = fc_vector_entrez, # Vector con Entrez IDs y LFC
                       circular = FALSE, colorEdge = TRUE, showCategory = 5,
                       node_label = "all", # Mostrar etiquetas de nodos (genes y vías)
                       cex_label_category = 1.2, cex_label_gene = 0.8,
                       title = paste0("Gene-Pathway Network (KEGG, ", comparison_prefix, ")")))
        dev.off()
      }, error = function(e) {
        message(paste("    ERROR al generar CnetPlot KEGG para", comparison_prefix, ": ", e$message))
        if(dev.cur() != 1) dev.off()
      })
    } else {
      message("  > No se pudo generar CnetPlot para KEGG (falta vector de LFC).")
    }
  } else {
    message(paste0("  > No se encontraron resultados KEGG para plotear para ", comparison_prefix, "."))
  }
  
  # DotPlot para Reactome
  if (!is.null(reactome_results) && nrow(as.data.frame(reactome_results)) > 0) {
    message("  Generando DotPlot para Reactome... (General)")
    png(paste0("DotPlot_Reactome_", comparison_prefix, ".png"), width = 900, height = 700, res = 120)
    print(dotplot(reactome_results, showCategory = 15, title = paste0("Reactome Pathway Enrichment (", comparison_prefix, ")")))
    dev.off()
    
    # Cnetplot para Reactome (requiere fc_vector_entrez)
    if (!is.null(fc_vector_entrez)) {
      message("  Generando CnetPlot para Reactome...")
      tryCatch({
        png(paste0("CnetPlot_Reactome_", comparison_prefix, ".png"), width = 1000, height = 900, res = 120)
        reactome_for_cnet <- reactome_results
        if (nrow(as.data.frame(reactome_results)) > 5) {
          reactome_for_cnet_df <- as.data.frame(reactome_results)
          reactome_for_cnet_df <- top_n(reactome_for_cnet_df, 5, wt = -p.adjust)
          reactome_for_cnet <- clusterProfiler::build_enrichResult(reactome_for_cnet_df)
        }
        print(cnetplot(reactome_results, # Usar el objeto enrichResult original
                       foldChange = fc_vector_entrez, # Vector con Entrez IDs y LFC
                       circular = FALSE, colorEdge = TRUE, showCategory = 5,
                       node_label = "all",
                       cex_label_category = 1.2, cex_label_gene = 0.8,
                       title = paste0("Gene-Pathway Network (Reactome, ", comparison_prefix, ")")))
        dev.off()
      }, error = function(e) {
        message(paste("    ERROR al generar CnetPlot Reactome para", comparison_prefix, ": ", e$message))
        if(dev.cur() != 1) dev.off()
      })
    } else {
      message("  > No se pudo generar CnetPlot para Reactome (falta vector de LFC).")
    }
  } else {
    message(paste0("  > No se encontraron resultados Reactome para plotear para ", comparison_prefix, "."))
  }
  
  # --- D) Pathway View (Solo para KEGG, opcional y más avanzado) ---
  # Esto requiere la instalación de 'pathview' y, a veces, Ghostscript.
  # Puede ser pesado de correr y generar muchos archivos si tienes muchas vías.
  # Se selecciona la vía KEGG más significativa para visualizar.
  if (!is.null(kegg_results) && nrow(as.data.frame(kegg_results)) > 0 && !is.null(fc_vector_entrez)) {
    message("  Generando Pathway View para la vía KEGG más enriquecida (si hay)...")
    # Ordenar por p.adjust y seleccionar la vía más significativa
    top_kegg_id_df <- as.data.frame(kegg_results) %>%
      arrange(p.adjust) %>%
      head(1)
    
    if (nrow(top_kegg_id_df) > 0) {
      top_kegg_id <- top_kegg_id_df$ID
      top_kegg_name <- top_kegg_id_df$Description
      
      message(paste0("  > Pathway View para KEGG ID: ", top_kegg_id, " (", top_kegg_name, ")"))
      tryCatch({
        # pathview guarda los archivos automáticamente en el directorio de trabajo
        pathview(gene.data = fc_vector_entrez,
                 pathway.id = top_kegg_id,
                 species = "hsa",
                 limit = list(gene=max(abs(fc_vector_entrez)), cpd=1), # Límites de color para genes y compuestos
                 out.suffix = paste0("KEGG_", comparison_prefix), # Añade sufijo para identificar el archivo
                 kegg.native = TRUE) # Usar la visualización nativa de KEGG (requiere internet)
        message(paste0("  > Pathway View para KEGG ID ", top_kegg_id, " guardado como 'hsa", top_kegg_id, ".", paste0("KEGG_", comparison_prefix), ".png'"))
      }, error = function(e) {
        message(paste("    ERROR al generar Pathway View para", comparison_prefix, " y KEGG ID ", top_kegg_id, ": ", e$message))
      })
    } else {
      message("  > No hay vías KEGG significativas para generar Pathway View.")
    }
  } else {
    message("  > No se puede generar Pathway View (KEGG) para esta comparación (faltan resultados KEGG o vector de LFC).")
  }
}

# --- EJECUTAR VISUALIZACIONES PARA CADA COMPARACIÓN ---
# Asegúrate de que los scripts de enriquecimiento (3, 5) y los de DEGs (2, 6)
# hayan generado los archivos necesarios con los nuevos nombres.
generate_enrichment_visualizations("WTC_vs_NTC")
generate_enrichment_visualizations("H1KOC_vs_NTC")
generate_enrichment_visualizations("H2KOC_vs_NTC")

message("Script 4_functional_visualization.R completado.")