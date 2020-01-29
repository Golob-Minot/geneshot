#!/bin/bash

# Run this script to test for invalid manifest file structures

# 0. Make sure you have a working installation of Nextflow and Docker
# 1. Download or clone this repository
# 2. Navigate to the repository folder
# 3. Run this script

for manifest_fp in data/mock.manifest.invalid1.csv data/mock.manifest.invalid2.csv; do

    echo $manifest_fp
    echo

    NXF_VER=19.10.0 nextflow run main.nf \
        -c nextflow.config \
        -profile testing \
        --manifest $manifest_fp \
        --output output/ \
        --hg_index data/hg_chr_21_bwa_index.tar.gz \
        -w work/ \
        -resume

done
