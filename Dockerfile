# Base image: R 4.4.2
FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive
ENV WORKDIR=/workspace

WORKDIR $WORKDIR

# System dependencies for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libpng-dev libjpeg-dev libtiff-dev libxt-dev \
    libblas-dev liblapack-dev libreadline-dev \
    libcairo2-dev libpango1.0-dev \
    bzip2 libbz2-dev xz-utils liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Install remotes
RUN Rscript -e 'install.packages("remotes", repos="https://cran.r-project.org")'

# CRAN packages with pinned versions
RUN Rscript -e 'remotes::install_version("lubridate", "1.9.4")' \
    && Rscript -e 'remotes::install_version("forcats", "1.0.0")' \
    && Rscript -e 'remotes::install_version("stringr", "1.5.1")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.1.4")' \
    && Rscript -e 'remotes::install_version("readr", "2.1.5")' \
    && Rscript -e 'remotes::install_version("tidyr", "1.3.1")' \
    && Rscript -e 'remotes::install_version("tibble", "3.2.1")' \
    && Rscript -e 'remotes::install_version("tidyverse", "2.0.0")' \
    && Rscript -e 'remotes::install_version("purrr", "1.0.4")' \
    && Rscript -e 'remotes::install_version("broom.mixed", "0.2.9.6")' \
    #&& Rscript -e 'remotes::install_version("broom", "1.0.8")' \
    && Rscript -e 'remotes::install_version("ggalluvial", "0.12.5")' \
    && Rscript -e 'remotes::install_version("rstatix", "0.7.2")' \
    && Rscript -e 'remotes::install_version("egg", "0.4.5")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("factoextra", "1.0.7")' \
    && Rscript -e 'remotes::install_version("circlize", "0.4.16")' \
    && Rscript -e 'remotes::install_version("RColorBrewer", "1.1-3")' \
    && Rscript -e 'remotes::install_version("tidytext", "0.4.2")' \
    && Rscript -e 'remotes::install_version("ggpubr", "0.6.0")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.5.2")'

# Install BiocManager for Bioconductor
RUN Rscript -e 'install.packages("BiocManager", repos="https://cran.r-project.org")' \
    && Rscript -e 'BiocManager::install(version="3.20", ask=FALSE)'

# Explicit Bioconductor packages pinned to sessionInfo versions
RUN Rscript -e 'BiocManager::install("CATALYST", version="1.30.2", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("HDCytoData", version="1.26.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("diffcyt", version="1.26.1", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("flowCore", version="2.18.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("ExperimentHub", version="2.14.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("AnnotationHub", version="3.14.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("BiocFileCache", version="2.14.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("SingleCellExperiment", version="1.28.1", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("SummarizedExperiment", version="1.36.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("Biobase", version="2.66.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("GenomicRanges", version="1.58.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("GenomeInfoDb", version="1.42.3", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("IRanges", version="2.40.1", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("S4Vectors", version="0.44.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("BiocGenerics", version="0.52.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("MatrixGenerics", version="1.18.1", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("matrixStats", version="1.5.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("BiocNeighbors", version="2.0.1", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("DelayedArray", version="0.32.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("XVector", version="0.46.0", ask=FALSE)' \
    && Rscript -e 'BiocManager::install("zlibbioc", version="1.52.0", ask=FALSE)'

# Copy project files
COPY . $WORKDIR

# Run the main script on container start
CMD ["Rscript", "source/final_scripts/run_all_scripts.R"]
