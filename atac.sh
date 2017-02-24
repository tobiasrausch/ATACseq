#!/bin/bash



SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Run analysis pipeline
${BASEDIR}/src/atac.sh $@
