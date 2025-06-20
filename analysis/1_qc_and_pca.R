# Cargar librerías
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr") # Para operaciones de datos
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer") # Para paletas de color en heatmaps
if (!requireNamespace("genefilter", quietly = TRUE)) BiocManager::install("genefilter")


library(DESeq2)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(genefilter) # Para rowVars

message("Todas las librerías cargadas.")

# Cargar datos de conteo
counts <- read.table("gene_counts.txt", header = TRUE, row.names = 1)
counts <- counts[, grep("Aligned", colnames(counts))]

# Definir los nombres de las muestras usando los nombres de los grupos y réplicas
# Asumiendo que el orden de las columnas en 'counts' es consistente (3 NTC, 3 WTC, 3 H1KOC, 3 H2KOC)
new_sample_names <- c(
  "NTC_R1", "NTC_R2", "NTC_R3",
  "WTC_R1", "WTC_R2", "WTC_R3",
  "H1KOC_R1", "H1KOC_R2", "H1KOC_R3",
  "H2KOC_R1", "H2KOC_R2", "H2KOC_R3"
)
colnames(counts) <- new_sample_names

# Definir condiciones
col_data <- data.frame(
  condition = factor(c(rep("NTC", 3), rep("WTC", 3), rep("H1KOC", 3), rep("H2KOC", 3)),
                     levels = c("NTC", "WTC", "H1KOC", "H2KOC")) # Usar los nuevos niveles
)
rownames(col_data) <- colnames(counts)

# Asegurarse de que el orden de las columnas en counts sea el mismo que el de las filas en col_data
counts <- counts[, rownames(col_data)]

message("Datos de conteo cargados y metadata preparada.")


# Crear objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col_data,
                              design = ~ condition)

# Pre-filtrado para reducir el tamaño del objeto y acelerar el cálculo
# Mantener solo genes con al menos 10 conteos en total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
message(paste0("Genes restantes después del pre-filtrado: ", nrow(dds)))

# Normalización de varianza (para PCA y Heatmaps)
message("Aplicando transformación de varianza estabilizada (vst)...")
vsd <- vst(dds, blind = TRUE) # blind = TRUE para análisis exploratorios

message("Objeto vst creado.")

# 1. PCA Plot
# Generar los datos para PCA
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Añadir los nombres de las réplicas como una columna en pcaData
pcaData$Replica <- rownames(pcaData)

# Plotear PCA
png("PCA_Plot_Samples_by_Condition.png", width = 800, height = 700, res = 120)
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = Replica), hjust = -0.1, vjust = -0.5, size = 3.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA de Muestras por Condición de Tratamiento") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
message("PCA generado y guardado en PCA_Plot_Samples_by_Condition.png")

# 2. Heatmap de los 1000 genes más variables
message("Generando heatmap de los 1000 genes más variables...")
# Calcular la varianza de los genes en los datos VST
# Se necesita genefilter::rowVars
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 1000)
heatmap_matrix <- assay(vsd)[topVarGenes, ]

# Ajustar la escala de la matriz para el heatmap (opcional, pheatmap lo hace por defecto si no se da breaks)
# pheatmap escala automáticamente si no se dan 'breaks'

png("heatmap_top_1000_variable_genes.png", width = 1000, height = 1200, res = 150)
pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  show_rownames = FALSE, # No mostrar todos los nombres de los genes
  cluster_cols = TRUE,
  show_colnames = TRUE, # Mostrar nombres de muestras (NTC_R1, WTC_R1, etc.)
  annotation_col = col_data,
  main = "Heatmap de los 1000 Genes Más Variables por Tratamiento",
  fontsize_col = 8,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Colores para expresión (ej. bajo, medio, alto)
  cutree_cols = 4 # Intenta cortar el dendrograma de columnas en 4 grupos si es deseado
)
dev.off()
message("Heatmap de los 1000 genes más variables generado y guardado en heatmap_top_1000_variable_genes.png")

message("Script 1_qc_and_pca.R completado.")