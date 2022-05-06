#!/bin/bash

set -euo pipefail

date
echo
echo "Running workflow '${WORKFLOW_REPO}' from ${PWD}"
echo

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${WORKFLOW_REPO}" \
    -params-file ._wb/tool/params.json \
    -resume

echo
date
echo Done
