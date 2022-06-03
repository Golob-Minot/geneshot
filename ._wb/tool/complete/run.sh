#!/bin/bash

set -euo pipefail

date
echo
echo "Running workflow '${TOOL_REPO}' from ${PWD}"
echo

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}" \
    -params-file ._wb/tool/params.json \
    -resume

echo
date
echo Done
