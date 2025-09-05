nextflow.enable.dsl=2

////////////////////////////////////////////////////
//                  INPUT PARAMETERS              //
////////////////////////////////////////////////////

params.samplesheet      = "/mnt/lustre/scratch/nlsas/home/ulc/co/sab/samplesheet.csv"
params.genome_dir       = "/mnt/lustre/scratch/nlsas/home/ulc/co/sab/genome"
params.genome_index     = "/mnt/lustre/scratch/nlsas/home/ulc/co/sab/genome_index/STAR_GRCh38_gencode_v48"
params.annotation       = "/mnt/lustre/scratch/nlsas/home/ulc/co/sab/genome/gencode.v48.primary_assembly.annotation.gtf"
params.output_dir       = "/mnt/lustre/scratch/nlsas/home/ulc/co/sab/output"
params.bam_dir          = "/mnt/lustre/scratch/nlsas/home/ulc/co/sab/bam_files"

////////////////////////////////////////////////////
//                CREATE DIRECTORIES              //
////////////////////////////////////////////////////

// Ensure output directories exist to avoid errors
[params.output_dir, params.bam_dir, params.genome_index].each { dir ->
    file(dir).mkdirs()
}

////////////////////////////////////////////////////
//                     PROCESSES                  //
////////////////////////////////////////////////////

// Quality trimming with fastp
process FASTP {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads1), path(reads2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trimmed.fq.gz"),
          path("${sample_id}_R2.trimmed.fq.gz"),
          path("${sample_id}.fastp.html"),
          path("${sample_id}.fastp.json")

    publishDir "${params.output_dir}/trimmed", mode: 'copy', pattern: "*_R*.trimmed.fq.gz"
    publishDir "${params.output_dir}/fastp_reports", mode: 'copy', pattern: "*.fastp.*"

    script:
    """
    fastp \
        --thread ${task.cpus} \
        --in1 $reads1 --in2 $reads2 \
        --out1 ${sample_id}_R1.trimmed.fq.gz \
        --out2 ${sample_id}_R2.trimmed.fq.gz \
        --html ${sample_id}.fastp.html \
        --json ${sample_id}.fastp.json
    """
}

// STAR genome index generation
process STAR_INDEX {
    tag "genome_index"

    input:
    path genome_fasta
    path annotation_gtf

    output:
    path "star_index"

    publishDir "${params.genome_index}", mode: 'copy'

    script:
    """
    mkdir -p star_index

    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles $genome_fasta \
         --sjdbGTFfile $annotation_gtf \
         --sjdbOverhang 100
    """
}

// STAR alignment
process STAR_ALIGN {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)
    path genome_index_dir

    output:
    tuple val(sample_id),
          path("${sample_id}.Aligned.sortedByCoord.out.bam"),
          path("${sample_id}.Log.final.out"),
          path("${sample_id}.Log.out"),
          path("${sample_id}.Log.progress.out"),
          path("${sample_id}.SJ.out.tab")

    publishDir params.bam_dir, mode: 'copy', pattern: "${sample_id}.*"

    script:
    """
    # Skip if BAM already exists
    if [ -f ${params.bam_dir}/${sample_id}.Aligned.sortedByCoord.out.bam ]; then
        echo "BAM for ${sample_id} already exists, skipping STAR_ALIGN"
        cp ${params.bam_dir}/${sample_id}.Aligned.sortedByCoord.out.bam .
        exit 0
    fi

    STAR --runThreadN ${task.cpus} \
         --genomeDir ${genome_index_dir} \
         --readFilesIn $r1 $r2 \
         --readFilesCommand zcat \
         --outFileNamePrefix "${sample_id}." \
         --outSAMtype BAM SortedByCoordinate

    # Check BAM creation
    if [ ! -s ${sample_id}.Aligned.sortedByCoord.out.bam ]; then
        echo "ERROR: BAM file not created properly" >&2
        exit 1
    fi
    """
}

// Gene counts using featureCounts
process FEATURECOUNTS {
    tag "all_samples"

    input:
    path bams
    path annotation_gtf

    output:
    path "gene_counts.txt"

    publishDir "${params.output_dir}/counts", mode: 'copy'

    script:
    """
    featureCounts -p --countReadPairs -T ${task.cpus} \
                  -a ${annotation_gtf} \
                  -o gene_counts.txt ${bams.join(" ")}
    """
}

////////////////////////////////////////////////////
//                     WORKFLOW                   //
////////////////////////////////////////////////////

workflow {

    // Read samplesheet and create sample tuples
    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.read1), file(row.read2)) }

    // Reference channels
    genome_fasta_ch   = Channel.fromPath("${params.genome_dir}/GRCh38.primary_assembly.genome.fa")
    annotation_gtf_ch = Channel.fromPath(params.annotation)

    // Reuse STAR index if exists
    def index_dir = file(params.genome_index)
    def hasIndex = index_dir.exists() && (index_dir.list() ?: []).size() > 5

    genome_index_ch = hasIndex \
        ? Channel.value(index_dir) \
        : STAR_INDEX(genome_fasta_ch, annotation_gtf_ch)

    // Trim FASTQ files
    fastp_out_ch = samples_ch | FASTP
    trimmed_ch   = fastp_out_ch.map { id, r1, r2, html, json -> tuple(id, r1, r2) }

    // Generate or reuse BAM files
    aligned_ch = trimmed_ch.map { id, r1, r2 ->
        def bam_file = file("${params.bam_dir}/${id}.Aligned.sortedByCoord.out.bam")
        if (bam_file.exists()) {
            tuple(id, bam_file)
        } else {
            STAR_ALIGN(tuple(id, r1, r2), genome_index_ch)
        }
    }

    // Extract BAMs and run featureCounts
    all_bams_ch = aligned_ch.map { it[1] }.collect()
    FEATURECOUNTS(all_bams_ch, annotation_gtf_ch)
}