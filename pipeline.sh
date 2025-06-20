#!/bin/bash

#SBATCH --job-name=pipeline_v3            # Nombre del trabajo
#SBATCH --output=salida_%j.log            # Archivo de salida (usa %j para ID de trabajo)
#SBATCH --error=error_%j.log              # Archivo de errores (usa %j para ID de trabajo)
#SBATCH --nodes=1                         # Numero de nodos
#SBATCH --ntasks=1                        # Numero de tareas
#SBATCH --cpus-per-task=8                 # Numero de CPUs por tarea
#SBATCH --mem=64G                         # Memoria por nodo (ajustar segun necesidad)
#SBATCH --time=12:00:00                    # Tiempo limite (HH:MM:SS)
#SBATCH --mail-type=ALL                   # Notificaciones por correo (END, FAIL, ALL)
#SBATCH --mail-user=santiago.ariza@udc.es

# Set up paths
GENOME_DIR="/mnt/lustre/scratch/nlsas/home/ulc/co/sab/genome"
GENOME_INDEX_DIR="/mnt/lustre/scratch/nlsas/home/ulc/co/sab/genome_index"
FASTQ_DIR="/mnt/lustre/scratch/nlsas/home/ulc/co/sab/HOSEPIC_EVs_RNA_seq"
OUTPUT_DIR="/mnt/lustre/scratch/nlsas/home/ulc/co/sab/output"

# Load required modules
module load star/2.7.9a
module load subread/2.0.5

# Build genome index
echo "Building genome index..."
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir "$GENOME_INDEX_DIR" \
     --genomeFastaFiles "$GENOME_DIR/GRCh38.primary_assembly.genome.fa" \
     --sjdbGTFfile "$GENOME_DIR/gencode.v48.primary_assembly.annotation.gtf" \
     --sjdbOverhang 100
echo "Genome index built."

# Align reads
echo "Starting read alignment..."
mkdir -p "$OUTPUT_DIR"

# Find all *_1.fq.gz files in subdirectories
find "$FASTQ_DIR" -type f -name "*_1.fq.gz" | while read -r SAMPLE_1; do
    # Derive the paired file name
    SAMPLE_2="${SAMPLE_1/_1.fq.gz/_2.fq.gz}"

    # Ensure the paired file exists
    if [[ -f "$SAMPLE_2" ]]; then
        # Extract basename without path or extensions
        BASENAME=$(basename "$SAMPLE_1" _1.fq.gz)
        echo "Processing pair: $SAMPLE_1 and $SAMPLE_2"

        STAR --runThreadN 8 \
             --genomeDir "$GENOME_INDEX_DIR" \
             --readFilesIn "$SAMPLE_1" "$SAMPLE_2" \
             --readFilesCommand zcat \
             --outFileNamePrefix "$OUTPUT_DIR/${BASENAME}." \
             --outSAMtype BAM SortedByCoordinate
    else
        echo "Warning: Paired file for $SAMPLE_1 not found. Skipping..."
    fi
done
echo "Read alignment completed."

# Count gene expression
echo "Counting gene expression..."

featureCounts -p \
              --countReadPairs \
              -T 8 \
              -a "$GENOME_DIR"/gencode.v48.primary_assembly.annotation.gtf \
              -o "$OUTPUT_DIR"/gene_counts.txt \
              "$OUTPUT_DIR"/*.Aligned.sortedByCoord.out.bam

echo "Gene counting completed."
echo "Ready for DEG analysis. Use DESeq2 in R."