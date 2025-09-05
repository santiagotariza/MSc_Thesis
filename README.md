# RNA-seq Analysis Pipeline

This repository contains a complete RNA-seq analysis pipeline, including **quality control, trimming, genome indexing, alignment, and gene quantification**.  
The pipeline is implemented in **Nextflow** (with an optional Bash script) for data processing, and in **R** for downstream exploratory analysis, normalization, differential expression, and functional interpretation.  

The workflow supports execution with **local HPC modules** and optionally uses **Singularity or Docker containers** for reproducibility.


## Features

- **Processing pipeline**:
  - Reads trimming with **fastp**
  - Alignment using **STAR**
  - Gene quantification with **featureCounts**
  - Optional reuse of existing STAR genome indexes
- **Analysis pipeline in R**:
  - Quality control and exploratory analysis
  - Normalization and batch-effect correction
  - Differential expression analysis
  - Functional enrichment: **ORA & GSEA**
  - Co-expression analysis with **lncRNAs**
  - ORF prediction and pathway integration
  - Automated reporting in RMarkdown
- **Reproducibility**:
  - Runs on HPC with modules or containers (Singularity/Docker)
  - Predefined folder structure for input and output


## Repository Structure

```
rnaseq-pipeline/
├── Dockerfile                  # Container definition (fastp, STAR, featureCounts, Nextflow, R/renv)
├── README.md                   # Project documentation
├── nextflow.config             # Nextflow configuration (modules + fallback containers)
├── main.nf                     # Nextflow workflow implementation
├── run_pipeline.sh             # Optional Bash pipeline
├── samplesheet.csv             # Example input samplesheet
├── genome/                     # Reference genome FASTA + annotation
│   ├── GRCh38.primary_assembly.genome.fa
│   └── gencode.v48.primary_assembly.annotation.gtf
├── genome_index/               # STAR genome index (generated/reused)
├── fastq/                      # Raw FASTQ files
├── fastp_reports/              # fastp QC reports (.html, .json)
├── output/                     # Outputs from the processing pipeline
│   ├── trimmed/                # Trimmed FASTQ files
│   ├── aligned/                # BAM files from STAR
│   └── counts/                 # Gene counts (featureCounts)
├── results/                    # Downstream R analysis results
└── scripts/                    # R analysis scripts
├── 00_run_all_analysis.R
├── 01_load_filter_data.R
├── 02_exploratory_qc_analysis.R
├── 03_normalize_transform.R
├── 04_exploratory_plots.R
├── 05_correct_batch_effect.R
├── 06_differential_expression.R
├── 07_ORA.R
├── 08_GSEA.R
├── 09_coexp_lncRNA_groups.R
├── 10_consolidate_coexp_results.R
├── 11_predict_lncRNA_ORFs.R
├── 12_integrate_orf_pathways.R
├── 13_prioritize_lncRNA_candidates.R
├── 14_render_report.R
└── rna_seq_report.Rmd
```


## Usage

### 1. Build Docker image

```bash
docker build -t rnaseq_pipeline:latest
```
On HPC systems without Docker, Singularity images will be used automatically if configured in nextflow.config.


### 2. Run Bash pipeline (stable option)

```bash
docker run --rm -v /path/to/project:/data rnaseq_pipeline:latest \
           /bin/bash /data/run_pipeline.sh
```

### 3. Run Nextflow pipeline (beta)

```bash
nextflow run main.nf -c nextflow.config
```

Or with Docker:

```bash
docker run --rm -v /path/to/project:/data rnaseq_pipeline:latest \
           /opt/nextflow/nextflow run /data/main.nf -c /data/nextflow.config
```


### 4. Run downstream R analysis

Once the gene_counts.txt file is generated (from featureCounts), place it in the project root directory. Then run:

```bash
Rscript scripts/00_run_all_analysis.R
```

This will sequentially execute all scripts in /scripts and generate results in /results, plus the final report rna_seq_report.html.

## Workflow Overview

```mermaid
flowchart TD

    subgraph Preprocessing["Pre-processing (Nextflow/Bash)"]
        A[FASTQ files] --> B[fastp : quality filtering & adapter trimming]
        B --> C[STAR : genome indexing & alignment]
        C --> D[featureCounts : generate raw gene counts]
        D --> E[gene_counts.txt]
    end

    subgraph Analysis["Downstream Analysis (R scripts)"]
        E --> F[01_load_filter_data.R : import count matrix, filter low-expressed genes]
        F --> G[02_exploratory_qc_analysis.R : assess library size, mapping rates, sample QC]
        G --> H[03_normalize_transform.R : apply normalization & variance-stabilizing transformations]
        H --> I[04_exploratory_plots.R : visualize sample clustering - PCA, heatmaps]
        I --> J[05_correct_batch_effect.R : adjust for technical covariates/batch effects]
        J --> K[06_differential_expression.R : identify differentially expressed genes]
        K --> L[07_ORA.R : test for enriched functional categories]
        K --> M[08_GSEA.R : run gene set enrichment with ranked lists]
        K --> N[09_coexp_lncRNA_groups.R : build lncRNA co-expression modules]
        N --> O[10_consolidate_coexp_results.R : summarize and merge module outputs]
        O --> P[11_predict_lncRNA_ORFs.R : detect potential coding ORFs within lncRNAs]
        P --> Q[12_integrate_orf_pathways.R : link ORFs to pathways & functional annotations]
        Q --> R[13_prioritize_lncRNA_candidates.R : rank lncRNAs by evidence & functional relevance]
        R --> S[14_render_report.R : compile results into final report]
        S --> T[rna_seq_report.html : interactive summary]
    end
```


## Configuration

- Nextflow parameters: nextflow.config
- Bash pipeline parameters: at the top of run_pipeline.sh
- R dependencies: managed with renv (renv.lock included in the repo)

⸻

## Notes

- Existing STAR genome indexes are automatically reused if detected.
- Containers ensure reproducibility, but the pipeline works with HPC modules as well.
- The R pipeline produces:
- QC plots, PCA, clustering
- Differential expression tables
- ORA/GSEA enrichment results
- Co-expression networks and ORF predictions
- A final RNA-seq analysis report (rna_seq_report.html)
