FROM rocker/r-ver:4.4.2
# ---- Base system ----
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# System dependencies for building R and most CRAN/Bioc packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gfortran wget curl ca-certificates gnupg \
    libreadline-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    libpng-dev libjpeg-dev libtiff-dev libblas-dev liblapack-dev \
    libxt-dev libx11-dev libfreetype6-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libhdf5-dev locales zlib1g-dev \
    && locale-gen en_US.UTF-8 \
    && rm -rf /var/lib/apt/lists/*

ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# ---- Install R 4.4.2 ----
RUN wget https://cran.r-project.org/src/base/R-4/R-4.4.2.tar.gz && \
    tar -xf R-4.4.2.tar.gz && \
    cd R-4.4.2 && \
    ./configure --enable-R-shlib --with-blas --with-lapack && \
    make -j$(nproc) && make install && \
    cd .. && rm -rf R-4.4.2*

# ---- Core R setup ----
RUN Rscript -e 'install.packages(c("remotes","BiocManager"), repos="https://cloud.r-project.org")'

# ---- Pin Bioconductor release ----
RUN Rscript -e 'BiocManager::install(version = "3.20", ask = FALSE)'

# ---- Install CRAN packages with versions ----
RUN Rscript -e 'remotes::install_version("lubridate", "1.9.4")' \
    && Rscript -e 'remotes::install_version("forcats", "1.0.0")' \
    && Rscript -e 'remotes::install_version("stringr", "1.5.1")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.1.4")' \
    && Rscript -e 'remotes::install_version("readr", "2.1.5")' \
    && Rscript -e 'remotes::install_version("tidyr", "1.3.1")' \
    && Rscript -e 'remotes::install_version("tibble", "3.2.1")' \
    && Rscript -e 'remotes::install_version("tidyverse", "2.0.0")' \
    && Rscript -e 'remotes::install_version("purrr", "1.0.4")' \
    && Rscript -e 'remotes::install_version("broom", "1.0.8")' \
    && Rscript -e 'remotes::install_version("broom.mixed", "0.2.9.6")' \
    && Rscript -e 'remotes::install_version("ggalluvial", "0.12.5")' \
    && Rscript -e 'remotes::install_version("rstatix", "0.7.2")' \
    && Rscript -e 'remotes::install_version("egg", "0.4.5")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("factoextra", "1.0.7")' \
    && Rscript -e 'remotes::install_version("RColorBrewer", "1.1-3")' \
    && Rscript -e 'remotes::install_version("tidytext", "0.4.2")' \
    && Rscript -e 'remotes::install_version("ggpubr", "0.6.0")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.5.2")' \
    && Rscript -e 'remotes::install_version("readxl", "1.4.5")' \
    && Rscript -e 'remotes::install_version("reshape2", "1.4.4")' \
    && Rscript -e 'remotes::install_version("glue", "1.8.0")'

# ---- Install Bioconductor packages ----
RUN Rscript -e 'BiocManager::install(c( \
    "ComplexHeatmap@2.22.0", \
    "CATALYST@1.30.2", \
    "HDCytoData@1.26.0", \
    "diffcyt@1.26.1", \
    "SingleCellExperiment@1.28.1", \
    "SummarizedExperiment@1.36.0", \
    "Biobase@2.66.0", \
    "GenomicRanges@1.58.0", \
    "IRanges@2.40.1", \
    "S4Vectors@0.44.0", \
    "MatrixGenerics@1.18.1", \
    "ExperimentHub@2.14.0", \
    "AnnotationHub@3.14.0", \
    "BiocFileCache@2.14.0", \
    "flowCore@2.18.0", \
    "cyCombine@0.2.19", \
    "iMUBAC@0.2.1", \
    "limma@3.62.2", \
    "FlowSOM@2.14.0", \
    "scran@1.34.0", \
    "scater@1.34.1" \
  ), ask=FALSE, update=FALSE)'

# ---- Optional: copy your project into the image ----
WORKDIR /workspace
COPY . /workspace

# ---- Default command ----
CMD ["R"]
