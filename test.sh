#!/bin/bash

set -e

# Run this script to test geneshot on the example data in the repository
# All execution will be done locally using Docker, so there is no need to
# set up HPC or cloud execution configurations

# 0. Make sure you have a working installation of Nextflow and Docker
# 1. Download or clone this repository
# 2. Navigate to the repository folder
# 3. Run this script

# Test with preprocessing and a formula
NXF_VER=19.10.0 nextflow run main.nf \
    -c nextflow.config \
    -profile testing \
    --manifest data/mock.manifest.csv \
    --preprocess_output output/preprocess_output \
    --output output1 \
    --hg_index data/hg_chr_21_bwa_index.tar.gz \
    --formula "label1 + label2" \
    --distance_threshold 0.1 \
    -w work/ \
    --noannot \
    --savereads \
    -resume

# Test with preprocessing and no formula
NXF_VER=19.10.0 nextflow run main.nf \
    -c nextflow.config \
    -profile testing \
    --manifest data/mock.manifest.csv \
    --preprocess_output output/preprocess_output \
    --output output2 \
    --hg_index data/hg_chr_21_bwa_index.tar.gz \
    --distance_threshold 0.1 \
    -w work/ \
    --noannot \
    -resume

# Test with formula and no preprocessing
NXF_VER=19.10.0 nextflow run main.nf \
    -c nextflow.config \
    -profile testing \
    --nopreprocess \
    --manifest data/mock.manifest.csv \
    --output output3 \
    --hg_index data/hg_chr_21_bwa_index.tar.gz \
    --formula "label1 + label2" \
    --distance_threshold 0.1 \
    -w work/ \
    --noannot \
    -resume

# Test with no formula and no preprocessing
NXF_VER=19.10.0 nextflow run main.nf \
    -c nextflow.config \
    -profile testing \
    --nopreprocess \
    --manifest data/mock.manifest.csv \
    --output output4 \
    --hg_index data/hg_chr_21_bwa_index.tar.gz \
    -w work/ \
    --noannot \
    -resume
