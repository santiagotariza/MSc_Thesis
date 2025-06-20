# Cargar librerías necesarias
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
if (!requireNamespace("apeglm", quietly = TRUE))
  BiocManager::install("apeglm")
if (!requireNamespace("EnhancedVolcano", quietly = TRUE))
  BiocManager::install("EnhancedVolcano")
if (!requireNamespace("ggplot2", quietly = TRUE)) # Para MA plot y Volcano
  install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) # Para manipulación de datos, aunque no estrictamente necesario aquí, es buena práctica
  install.packages("dplyr")
if (!requireNamespace("tibble", quietly = TRUE)) # Para rownames_to_column
  install.packages("tibble")


library(DESeq2)
library(apeglm)
library(EnhancedVolcano)
library(ggplot2)
library(dplyr)
library(tibble) # Asegurarse de que esté cargado para rownames_to_column

message("Librerías cargadas. Iniciando análisis de expresión diferencial para GSEA...")

# Cargar datos de conteo (asegurarse de que estén preprocesados)
counts <- read.table("gene_counts.txt", header = TRUE, row.names = 1)

# Filtrar columnas relevantes con "Aligned" en sus nombres
counts <- counts[, grep("Aligned", colnames(counts))]

# Definir los nombres de las muestras usando los nombres de los grupos y réplicas
new_sample_names <- c(
  "NTC_R1", "NTC_R2", "NTC_R3",
  "WTC_R1", "WTC_R2", "WTC_R3",
  "H1KOC_R1", "H1KOC_R2", "H1KOC_R3",
  "H2KOC_R1", "H2KOC_R2", "H2KOC_R3"
)
colnames(counts) <- new_sample_names

# Definir condiciones (usando los nuevos nombres de grupos)
col_data <- data.frame(
  condition = factor(c(rep("NTC", 3), rep("WTC", 3), rep("H1KOC", 3), rep("H2KOC", 3)),
                     levels = c("NTC", "WTC", "H1KOC", "H2KOC"))
)
rownames(col_data) <- colnames(counts)

# Asegurarse de que el orden de las columnas en counts sea el mismo que el de las filas en col_data
counts <- counts[, rownames(col_data)]

message("Datos de conteo cargados y metadata preparada.")


# Crear objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col_data,
                              design = ~ condition)

# Pre-filtrado (opcional pero recomendado)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
message(paste0("Genes restantes después del pre-filtrado: ", nrow(dds)))

# Ejecutar el pipeline DESeq2
message("Corriendo DESeq2...")
dds <- DESeq(dds)
message("DESeq2 completado.")


# Función para realizar la comparación, encogimiento de LFC y guardar resultados completos
analyze_and_save_degs <- function(dds_obj, contrast_levels, output_prefix) {
  message(paste0("  Analizando contraste para GSEA: ", paste(contrast_levels, collapse = " vs ")))
  
  # Obtener resultados sin encogimiento (se usa para el LFC Shrink)
  res <- results(dds_obj, contrast = c("condition", contrast_levels[1], contrast_levels[2]), alpha = 0.05)
  
  # Aplicar encogimiento de LFC (importante para GSEA)
  res_lfc_shrink <- lfcShrink(dds_obj, contrast = c("condition", contrast_levels[1], contrast_levels[2]), res = res, type = "apeglm")
  
  # Convertir a data.frame y añadir columna de GeneID
  res_df <- as.data.frame(res_lfc_shrink) %>%
    rownames_to_column(var = "GeneID")
  
  # Guardar los resultados completos
  write.csv(res_df, paste0("DESeq2_results_full_", output_prefix, ".csv"), row.names = FALSE)
  message(paste0("  > Resultados completos guardados en DESeq2_results_full_", output_prefix, ".csv"))
  
  # Opcional: Generar Volcano plot (solo para verificación visual)
  tryCatch({
    png(paste0("Volcano_plot_GSEA_prep_", output_prefix, ".png"), width = 900, height = 800, res = 120)
    print(EnhancedVolcano(res_lfc_shrink,
                          lab = rownames(res_lfc_shrink),
                          x = 'log2FoldChange',
                          y = 'padj',
                          title = paste("Volcano Plot (GSEA Prep):", output_prefix),
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
    message(paste("    Volcano plot (GSEA prep) guardado en:", paste0("Volcano_plot_GSEA_prep_", output_prefix, ".png")))
  }, error = function(e) {
    message(paste("    ERROR al generar Volcano plot (GSEA prep) para", output_prefix, ": ", e$message))
    if(dev.cur() != 1) dev.off() # Asegurarse de cerrar el dispositivo gráfico
  })
  
  return(res_lfc_shrink) # Retorna los resultados completos
}

# Realizar análisis para cada comparación (¡AJUSTADAS LAS COMPARACIONES AQUÍ!)
res_WTC_vs_NTC_full <- analyze_and_save_degs(dds, c("WTC", "NTC"), "WTC_vs_NTC")
res_H1KOC_vs_NTC_full <- analyze_and_save_degs(dds, c("H1KOC", "NTC"), "H1KOC_vs_NTC")
res_H2KOC_vs_NTC_full <- analyze_and_save_degs(dds, c("H2KOC", "NTC"), "H2KOC_vs_NTC")

message("Análisis de expresión diferencial para GSEA completado. Archivos CSV y plots generados.")

message("Script 6_differential_expression_for_GSEA.R completado.")