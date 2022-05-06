#!/bin/bash

set -em

echo "Setting up nextflow.config"

echo """
docker.enabled = true
report.enabled = true
trace.enabled = true
""" > nextflow.config

# If the Tower token was provided
if [ -z ${TOWER_ACCESS_TOKEN} ]; then
    echo Tower token is not set
else
    echo """
tower {
  accessToken = '${TOWER_ACCESS_TOKEN}'
  enabled = true
}
""" >> nextflow.config
fi

cat nextflow.config
echo

# Disable ANSI logging
export NXF_ANSI_LOG=false

# Print the Nextflow version being used
echo "Nextflow Version: ${NXF_VER}"
echo

# Execute the tool in the local environment
echo "Starting tool"
echo

# Start the tool
/bin/bash ._wb/helpers/run_tool
