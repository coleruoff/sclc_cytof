FROM rocker/r-ver:4.4.2

LABEL maintainer="Cole Ruoff <coleruoff16@gmail.com>"
LABEL description="Reproducible R environment for my SCLC CTCs CyTOF project"

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    make \
    gcc \
    g++ \
    gfortran \
    libpng-dev \
    libjpeg-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install renv (needed to restore packages)
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"

# Set working directory inside container
WORKDIR /sclc_cytof

# Copy your project into the container
COPY . /sclc_cytof

# Restore all R packages exactly as recorded in renv.lock
RUN R -e "renv::restore()"

# Default command: run R scripts
ENTRYPOINT ["Rscript"]