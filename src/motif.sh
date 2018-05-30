#!/bin/bash

if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: $0 <peaks.tsv> <output prefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Activate environment
export PATH=${BASEDIR}/../bin/bin:${PATH}
source activate ${BASEDIR}/../bin/envs/atac2

# CMD parameters
PEAKS=${1}
OUTP=${2}

# Run motif finding
mkdir ${OUTP}_motifs
findMotifsGenome.pl ${PEAKS} hg19 ${OUTP}_motifs/ -size 50 -mask

# Deactivate environment
source deactivate
