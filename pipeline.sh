#!/bin/bash

#SBATCH --job-name=RNAseq_FASTP
#SBATCH --output=stdout_%j.log
#SBATCH --error=stderr_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=santiago.ariza@udc.es

# Paths and parameters
GENOME_DIR="/mnt/lustre/scratch/nlsas/home/ulc/co/sab/genome"
GENOME_INDEX_DIR="/mnt/lustre/scratch/nlsas/home/ulc/co/sab/genome_index"
FASTQ_DIR="/mnt/lustre/scratch/nlsas/home/ulc/co/sab/HOSEPIC_EVs_RNA_seq"
OUTPUT_DIR="/mnt/lustre/scratch/nlsas/home/ulc/co/sab/output"

TRIMMED_DIR="$OUTPUT_DIR/trimmed"
FASTP_HTML="$OUTPUT_DIR/fastp_reports"
ALIGN_DIR="$OUTPUT_DIR/aligned"
COUNT_DIR="$OUTPUT_DIR/counts"

ERROR_LOG="$OUTPUT_DIR/pipeline_errors.log"
THREADS=8

# Load modules
module load fastp/0.22.0
module load star/2.7.9a
module load subread/2.0.5

# Create output directories
mkdir -p "$TRIMMED_DIR" "$FASTP_HTML" "$ALIGN_DIR" "$COUNT_DIR"
echo "" > "$ERROR_LOG"

# Step 1: fastp (QC + trimming)
echo "Step 1: Running fastp for trimming and quality filtering..."
find "$FASTQ_DIR" -type f -name "*_1.fq.gz" | while read -r SAMPLE_1; do
    SAMPLE_2="${SAMPLE_1/_1.fq.gz/_2.fq.gz}"

    if [[ -f "$SAMPLE_2" ]]; then
        BASENAME=$(basename "$SAMPLE_1" _1.fq.gz)

        OUT_R1="$TRIMMED_DIR/${BASENAME}_R1.trimmed.fq.gz"
        OUT_R2="$TRIMMED_DIR/${BASENAME}_R2.trimmed.fq.gz"
        HTML_OUT="$FASTP_HTML/${BASENAME}.fastp.html"
        JSON_OUT="$FASTP_HTML/${BASENAME}.fastp.json"

        fastp \
            --thread $THREADS \
            --in1 "$SAMPLE_1" --in2 "$SAMPLE_2" \
            --out1 "$OUT_R1" --out2 "$OUT_R2" \
            --detect_adapter_for_pe \
            --qualified_quality_phred 15 \
            --length_required 36 \
            --html "$HTML_OUT" \
            --json "$JSON_OUT"

        if [[ $? -ne 0 ]]; then
            echo "[ERROR] fastp failed for sample: $BASENAME" >> "$ERROR_LOG"
        fi
    else
        echo "[WARNING] Paired file missing for: $SAMPLE_1" >> "$ERROR_LOG"
    fi
done

# Step 2: Genome indexing (only if needed)
if [ ! -f "$GENOME_INDEX_DIR/SA" ]; then
    echo "Step 2: Building STAR genome index..."
    STAR --runThreadN $THREADS \
         --runMode genomeGenerate \
         --genomeDir "$GENOME_INDEX_DIR" \
         --genomeFastaFiles "$GENOME_DIR/GRCh38.primary_assembly.genome.fa" \
         --sjdbGTFfile "$GENOME_DIR/gencode.v48.primary_assembly.annotation.gtf" \
         --sjdbOverhang 100

    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Genome indexing failed" >> "$ERROR_LOG"
        exit 1
    fi
else
    echo "Step 2: Genome index already exists. Skipping indexing."
fi

# Step 3: STAR alignment
echo "Step 3: Aligning reads with STAR..."
find "$TRIMMED_DIR" -name "*_R1.trimmed.fq.gz" | while read -r TRIM_1; do
    TRIM_2="${TRIM_1/_R1.trimmed.fq.gz/_R2.trimmed.fq.gz}"

    if [[ -f "$TRIM_2" ]]; then
        BASENAME=$(basename "$TRIM_1" _R1.trimmed.fq.gz)

        STAR --runThreadN $THREADS \
             --genomeDir "$GENOME_INDEX_DIR" \
             --readFilesIn "$TRIM_1" "$TRIM_2" \
             --readFilesCommand zcat \
             --outFileNamePrefix "$ALIGN_DIR/${BASENAME}." \
             --outSAMtype BAM SortedByCoordinate

        if [[ $? -ne 0 ]]; then
            echo "[ERROR] STAR alignment failed for sample: $BASENAME" >> "$ERROR_LOG"
        fi
    else
        echo "[WARNING] Paired trimmed file missing for: $TRIM_1" >> "$ERROR_LOG"
    fi
done

# Step 4: Gene quantification with featureCounts
echo "Step 4: Counting gene expression with featureCounts..."
featureCounts -p --countReadPairs -T $THREADS \
              -a "$GENOME_DIR/gencode.v48.primary_assembly.annotation.gtf" \
              -o "$COUNT_DIR/gene_counts.txt" \
              "$ALIGN_DIR"/*.Aligned.sortedByCoord.out.bam

if [[ $? -ne 0 ]]; then
    echo "[ERROR] featureCounts failed" >> "$ERROR_LOG"
    exit 1
fi

# Final check
if [ -s "$ERROR_LOG" ]; then
    echo "Pipeline completed with warnings/errors. See $ERROR_LOG for details."
else
    echo "Pipeline completed successfully with no errors."
fi