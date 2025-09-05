# ========================================
# Combined RNA-seq Analysis & Processing Container
# Includes R 4.5.1 with renv, fastp, STAR, featureCounts, and Nextflow
# ========================================

# Base image
FROM ubuntu:22.04

# -----------------------------
# Environment variables
# -----------------------------
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/nextflow:${PATH}"

# -----------------------------
# Install system dependencies
# -----------------------------
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    git \
    gzip \
    zlib1g-dev \
    build-essential \
    openjdk-11-jre-headless \
    python3 \
    python3-pip \
    locales \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglpk-dev \
    libxt-dev \
    libcairo2-dev \
    libx11-dev \
    && locale-gen en_US.UTF-8 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------
# Install R 4.5.1
# -----------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base=4.5.1* \
    r-base-dev=4.5.1* \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------
# Install renv for package management
# -----------------------------
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"



# --- If your system already has the next packages then you may skip the following installation steps ---


# -----------------------------
# Install fastp
# -----------------------------
RUN wget -O /usr/local/bin/fastp https://github.com/OpenGene/fastp/releases/download/v0.22.0/fastp \
    && chmod +x /usr/local/bin/fastp

# -----------------------------
# Install STAR
# -----------------------------
RUN wget -O /tmp/STAR.tar.gz https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz \
    && tar -xzf /tmp/STAR.tar.gz -C /tmp \
    && cd /tmp/STAR-2.7.9a/source \
    && make STAR \
    && cp STAR /usr/local/bin/ \
    && cd / \
    && rm -rf /tmp/STAR*

# -----------------------------
# Install featureCounts (Subread)
# -----------------------------
RUN wget -O /tmp/subread-2.0.5-Linux-x86_64.tar.gz https://sourceforge.net/projects/subread/files/subread-2.0.5/subread-2.0.5-Linux-x86_64.tar.gz/download \
    && tar -xzf /tmp/subread-2.0.5-Linux-x86_64.tar.gz -C /opt/ \
    && ln -s /opt/subread-2.0.5-Linux-x86_64/bin/featureCounts /usr/local/bin/featureCounts \
    && rm -rf /tmp/subread*

# -----------------------------
# Install Nextflow (Optional: only use it if you plan to run the Nextflow pipeline version)
# -----------------------------
#RUN curl -s https://get.nextflow.io | bash \
#    && mv nextflow /opt/nextflow/nextflow \
#    && chmod +x /opt/nextflow/nextflow

# -----------------------------
# Set working directories
# -----------------------------
RUN mkdir -p /data/scripts /data/results /data/fastp_reports /data/genome /data/genome_index /data/fastq

# Set working directory
WORKDIR /data

# Copy project files (including renv.lock if available)
COPY . .

# Restore R environment if renv.lock is present
RUN R -e "if (file.exists('renv.lock')) renv::restore(confirm = FALSE)"

# Default command
CMD ["R"]