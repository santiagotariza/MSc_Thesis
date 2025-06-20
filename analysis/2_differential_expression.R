# Instalar librerías si no están instaladas
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("apeglm", quietly = TRUE))
  BiocManager::install("apeglm")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
if (!requireNamespace("EnhancedVolcano", quietly = TRUE))
  BiocManager::install("EnhancedVolcano")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("tibble", quietly = TRUE)) # Para rownames_to_column
  install.packages("tibble")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) # Para mapeo de IDs
  BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) # Para mapIds
  BiocManager::install("AnnotationDbi")


# Cargar librerías
library(DESeq2)
library(apeglm) # Para lfcShrink
library(EnhancedVolcano) # Para Volcano plots avanzados
library(ggplot2) # Para MA plot y Volcano
library(dplyr) # Para manipulación de datos
library(tibble) # Para rownames_to_column
library(org.Hs.eg.db) # Base de datos de anotación para Homo sapiens
library(AnnotationDbi) # Para la función mapIds

message("Librerías cargadas. Iniciando análisis de expresión diferencial...")

# Cargar datos de conteo
counts <- read.table("gene_counts.txt", header = TRUE, row.names = 1)

# Filtrar columnas relevantes con "Aligned" en sus nombres
counts <- counts[, grep("Aligned", colnames(counts))]

# Simplificar nombres de las muestras y asegurar el orden correcto
# Basado en la confirmación:
# NTC -> G1 (No-Treatment Control)
# WTC -> G2 (Wild-Type Extracellular Vesicles Treated)
# H1KOC -> G3 (HMGB1 Knockout Extracellular Vesicles Treated)
# H2KOC -> G4 (HMGB2 Knockout Extracellular Vesicles Treated)

# Renombrar las columnas para que coincidan con los grupos confirmados
# Asumiendo que el orden de las columnas es NTC1, NTC2, NTC3, WTC1, WTC2, WTC3, H1KOC1, H1KOC2, H1KOC3, H2KOC1, H2KOC2, H2KOC3
new_colnames <- c("NTC1", "NTC2", "NTC3", "WTC1", "WTC2", "WTC3", "H1KOC1", "H1KOC2", "H1KOC3", "H2KOC1", "H2KOC2", "H2KOC3")
colnames(counts) <- new_colnames

# Definir la información de las columnas (condiciones de las muestras)
col_data <- data.frame(
  condition = factor(c(rep("NTC", 3), rep("WTC", 3), rep("H1KOC", 3), rep("H2KOC", 3)),
                     levels = c("NTC", "WTC", "H1KOC", "H2KOC")) # Define el orden y NTC como referencia
)
rownames(col_data) <- colnames(counts) # Asegúrate de que los nombres de las filas coincidan con los nombres de las columnas de counts

# Asegurarse de que el orden de las columnas en counts sea el mismo que el de las filas en col_data
counts <- counts[, rownames(col_data)]

# Filtrar genes con conteos muy bajos
# Se eliminan genes con menos de 10 conteos en al menos 3 muestras (número mínimo de muestras por grupo)
keep <- rowSums(counts >= 10) >= 3 # Ajusta el '3' según el tamaño mínimo de tu grupo más pequeño
counts_filtered <- counts[keep, ]
message(paste0("Genes antes del filtrado: ", nrow(counts), ". Genes después del filtrado: ", nrow(counts_filtered), "."))

# Crear objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = col_data,
                              design = ~ condition)

# Ejecutar el análisis DESeq2
message("Corriendo el análisis DESeq2. Esto puede tardar unos minutos...")
dds <- DESeq(dds)
message("Análisis DESeq2 completado.")

# --- Función para realizar el análisis de expresión diferencial y guardar los resultados ---
analyze_comparison <- function(dds_obj, contrast_levels, output_prefix) {
  message(paste0("  Analizando contraste: ", paste(contrast_levels, collapse = " vs ")))
  
  # 1. Realizar el contraste inicial con results()
  # 'contrast_levels[1]' es el numerador (ej. WTC), 'contrast_levels[2]' es el denominador (ej. NTC)
  res <- results(dds_obj, contrast = c("condition", contrast_levels[1], contrast_levels[2]), alpha = 0.05)
  
  # 2. Preparar el nombre del coeficiente para lfcShrink con apeglm
  # DESeq2 crea coeficientes con nombres como "condition_FACTOR_vs_REFERENCE"
  coef_name_expected <- paste0("condition_", contrast_levels[1], "_vs_", contrast_levels[2])
  
  # Verificar si el coeficiente existe antes de usarlo
  all_coefs <- resultsNames(dds_obj) # Obtener todos los nombres de coeficientes del modelo
  
  res_lfc_shrink <- NULL # Inicializar
  if (coef_name_expected %in% all_coefs) {
    message(paste0("    Aplicando shrinkage con apeglm para el coeficiente: '", coef_name_expected, "'"))
    res_lfc_shrink <- lfcShrink(dds_obj, coef = coef_name_expected, type = "apeglm")
  } else {
    message(paste0("    ADVERTENCIA: Coeficiente esperado '", coef_name_expected, "' no encontrado en resultsNames(dds_obj)."))
    message("    Coeficientes disponibles son: ", paste(all_coefs, collapse = ", "))
    message("    No se puede aplicar apeglm shrinkage sin un nombre de coeficiente válido.")
    message("    Los resultados se guardarán sin shrinkage de LFC.")
    res_lfc_shrink <- res # Usar los resultados no encogidos si apeglm no puede aplicarse
  }
  
  # Convertir a un data frame y añadir GeneID
  res_df <- as.data.frame(res_lfc_shrink) %>%
    rownames_to_column(var = "GeneID")
  
  # Mapear IDs de Ensembl a Gene Symbols
  # Limpiar IDs de Ensembl de la versión (e.g., ENSG0000012345.10 -> ENSG0000012345)
  res_df$ensembl_id_clean <- sub("\\..*", "", res_df$GeneID)
  
  # Mapear a Gene Symbols usando AnnotationDbi::mapIds
  suppressMessages({ # Para evitar mensajes sobre multiVals
    res_df$GeneSymbol <- mapIds(org.Hs.eg.db,
                                keys = res_df$ensembl_id_clean,
                                column = "SYMBOL",
                                keytype = "ENSEMBL",
                                multiVals = "first") # Tomar el primer símbolo si hay múltiples mapeos
  })
  
  # Guardar los resultados completos antes de filtrar DEGs
  output_full_path <- paste0("DESeq2_results_full_", output_prefix, ".csv")
  write.csv(res_df, output_full_path, row.names = FALSE)
  message(paste0("    Resultados DESeq2 completos guardados en: ", output_full_path))
  
  # Filtrar DEGs: padj < 0.05 y |log2FoldChange| > 1
  degs <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
  
  # Ordenar por p.adjust
  degs <- degs[order(degs$padj), ]
  
  # Guardar los DEGs
  output_degs_path <- paste0("DEGs_", output_prefix, ".csv")
  if (nrow(degs) > 0) {
    write.csv(degs, output_degs_path, row.names = FALSE)
    message(paste0("    DEGs guardados en: ", output_degs_path, " (", nrow(degs), " genes)"))
  } else {
    message(paste0("    No se encontraron DEGs significativos para ", output_prefix, "."))
  }
  
  # Generar MA Plot
  ma_plot_path <- paste0("MA_plot_", output_prefix, ".png")
  png(ma_plot_path, width = 800, height = 700, res = 120)
  plotMA(res_lfc_shrink, main = paste0("MA Plot: ", paste(contrast_levels, collapse = " vs ")), ylim = c(-5, 5))
  dev.off()
  message(paste0("    MA plot guardado en: ", ma_plot_path))
  
  # Generar Volcano Plot usando EnhancedVolcano
  # AHORA USAMOS dplyr::select EXPLÍCITAMENTE
  volcano_plot_data <- as.data.frame(res_lfc_shrink) %>%
    rownames_to_column(var = "GeneID") %>%
    left_join(res_df %>% dplyr::select(GeneID, GeneSymbol), by = "GeneID") %>% # CORRECCIÓN AQUÍ
    column_to_rownames(var = "GeneID")
  
  # Asegurarse de que las etiquetas de genes no sean NA
  volcano_plot_data$GeneSymbol[is.na(volcano_plot_data$GeneSymbol)] <- rownames(volcano_plot_data)[is.na(volcano_plot_data$GeneSymbol)]
  
  volcano_plot_path <- paste0("Volcano_plot_", output_prefix, ".png")
  png(volcano_plot_path, width = 1000, height = 1000, res = 150)
  tryCatch({
    print(EnhancedVolcano(volcano_plot_data,
                          lab = volcano_plot_data$GeneSymbol, # Usar GeneSymbols para las etiquetas
                          x = 'log2FoldChange',
                          y = 'padj',
                          title = paste("Volcano Plot:", paste(contrast_levels, collapse = " vs ")),
                          pCutoff = 0.05,
                          FCcutoff = 1.0,
                          pointSize = 1.5,
                          labSize = 3,
                          legendLabSize = 12,
                          legendIconSize = 4.0,
                          drawConnectors = TRUE,
                          widthConnectors = 0.5,
                          colConnectors = 'grey30',
                          subtitle = NULL)) 
    dev.off()
    message(paste("    Volcano plot guardado en:", volcano_plot_path))
  }, error = function(e) {
    message(paste("    ERROR al generar Volcano plot para", output_prefix, ": ", e$message))
    if(dev.cur() != 1) dev.off() 
  })
  
  return(res_lfc_shrink) # Retorna los resultados encogidos
}

# --- Realizar análisis para cada comparación ---
# ¡AJUSTADAS LAS COMPARACIONES A LA NUEVA NOMENCLATURA!
message("\nIniciando análisis para las comparaciones específicas...")
res_WTC_vs_NTC_full <- analyze_comparison(dds, c("WTC", "NTC"), "WTC_vs_NTC")
res_H1KOC_vs_NTC_full <- analyze_comparison(dds, c("H1KOC", "NTC"), "H1KOC_vs_NTC")
res_H2KOC_vs_NTC_full <- analyze_comparison(dds, c("H2KOC", "NTC"), "H2KOC_vs_NTC")

message("\nTodos los análisis de expresión diferencial completados. Archivos CSV y plots generados.")