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

# Test with the gene catalog made in a previous round
NXF_VER=19.10.0 nextflow run main.nf \
    -c nextflow.config \
    -profile testing \
    --gene_fasta output1/ref/genes.fasta.gz \
    --nopreprocess \
    --manifest data/mock.manifest.csv \
    --output output5 \
    --hg_index data/hg_chr_21_bwa_index.tar.gz \
    -w work/ \
    --noannot \
    -resume


# Test with the gene catalog made in a previous round and whole genome alignment
NXF_VER=19.10.0 nextflow run main.nf \
    -c nextflow.config \
    -profile testing \
    --gene_fasta output1/ref/genes.fasta.gz \
    --nopreprocess \
    --manifest data/mock.manifest.csv \
    --output output6 \
    --hg_index data/hg_chr_21_bwa_index.tar.gz \
    -w work/ \
    --eggnog_db false \
    --eggnog_dmnd false \
    --taxonomic_dmnd false \
    --ref_genome_fasta data/reference_genomes.fasta.gz \
    --ref_genome_csv data/reference_genomes.csv.gz \
    -resume

# Test with de novo assembly and whole genome alignment
NXF_VER=19.10.0 nextflow run main.nf \
    -c nextflow.config \
    -profile testing \
    --nopreprocess \
    --manifest data/mock.manifest.csv \
    --output output7 \
    --hg_index data/hg_chr_21_bwa_index.tar.gz \
    --distance_threshold 0.1 \
    -w work/ \
    --eggnog_db false \
    --eggnog_dmnd false \
    --taxonomic_dmnd false \
    --ref_genome_fasta data/reference_genomes.fasta.gz \
    --ref_genome_csv data/reference_genomes.csv.gz \
    -resume
