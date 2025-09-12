# ========================================
# RNA-seq Analysis & Processing Container
# Includes R 4.5.1 with renv, fastp, STAR, featureCounts, and optionally Nextflow
# ========================================

# Base image
FROM ubuntu:22.04

# -----------------------------
# Environment setup
# -----------------------------
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=en_US.UTF-8
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
    libncurses5-dev \
    libbz2-dev \
    xz-utils \
    liblzma-dev \
    && locale-gen en_US.UTF-8 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# -----------------------------
# Install R (latest stable from CRAN)
# -----------------------------
RUN apt-get update && apt-get install -y --no-install-recommends software-properties-common dirmngr gnupg apt-transport-https ca-certificates \
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-base r-base-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------
# Install renv for R package management
# -----------------------------
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"

# -----------------------------
# Install fastp
# -----------------------------
RUN wget -O /usr/local/bin/fastp http://opengene.org/fastp/fastp.0.22.0 \
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
# Install featureCounts (from Subread)
# -----------------------------
RUN wget -O /tmp/subread.tar.gz https://sourceforge.net/projects/subread/files/subread-2.0.5/subread-2.0.5-Linux-x86_64.tar.gz/download \
    && tar -xzf /tmp/subread.tar.gz -C /opt/ \
    && ln -s /opt/subread-2.0.5-Linux-x86_64/bin/featureCounts /usr/local/bin/featureCounts \
    && rm -rf /tmp/subread*

# -----------------------------
# Install Nextflow (optional)
# -----------------------------
# RUN curl -s https://get.nextflow.io | bash \
#     && mkdir -p /opt/nextflow \
#     && mv nextflow /opt/nextflow/nextflow \
#     && chmod +x /opt/nextflow/nextflow

# -----------------------------
# Set working directory in container
# -----------------------------
WORKDIR /app

# Copiar archivos del proyecto (incluyendo scripts/ y renv.lock)
COPY . .

# Inicializar y restaurar entorno renv
RUN R -e "renv::init(bare = TRUE); renv::restore(confirm = FALSE)"

# -----------------------------
# Default command
# -----------------------------
CMD ["R"]