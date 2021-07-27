container__find_cags = "quay.io/fhcrc-microbiome/python-pandas:sklearn"

// Default options
params.distance_threshold = 0.15
params.distance_metric = "cosine"
params.min_contig_size = 3
params.min_contig_depth = 5

// Group together genes by co-abundance
process makeCAGs {
    tag "Group genes by co-abundance"
    container "${container__find_cags}"
    label "mem_veryhigh"

    input:
    path assembly_hdf
    path alignment_hdf

    output:
    file "CAGs.csv.gz"
    file "CAGs.time.json"

    """#!/bin/bash

set -Eeuo pipefail

export NUMEXPR_MAX_THREADS=${task.cpus}

make_cags.py \
    --assembly-hdf ${assembly_hdf} \
    --abundance-hdf ${alignment_hdf} \
    --metric ${params.distance_metric} \
    --threshold ${params.distance_threshold} \
    --min-contig-size ${params.min_contig_size} \
    --min-contig-depth ${params.min_contig_depth} \
    --min-specimens ${params.min_specimens} \
    --output-folder ./ \
    --output-prefix CAGs \
    --processes ${task.cpus}

    """
}
