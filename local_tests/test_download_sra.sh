#!/bin/bash

# Run this script to test the `download_sra.nf` workflow in the `geneshot` repository
# All execution will be done locally using Docker, so there is no need to
# set up HPC or cloud execution configurations

# 0. Make sure you have a working installation of Nextflow and Docker
# 1. Download or clone this repository
# 2. Navigate to the repository folder
# 3. Run this script

# Test with a whole BioProject of paired-end reads
NXF_VER=20.01.0 nextflow run download_sra.nf \
    -c nextflow.config.sample \
    -profile testing \
    -with-docker ubuntu:20.04 \
    --accession PRJNA541981 \
    --output output_download_1 \
    --metadata \
    -w work/ \
    -resume
