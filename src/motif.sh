#!/bin/bash

if [ $# -ne 3 ]
then
    echo ""
    echo "Usage: $0 <hg19|mm10> <peaks.tsv> <output prefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Activate environment
export PATH=${BASEDIR}/../bin/bin:${PATH}
source activate ${BASEDIR}/../bin/envs/atac2

# CMD parameters
ATYPE=${1}
PEAKS=${2}
OUTP=${3}

# Run motif finding
mkdir ${OUTP}_motifs
findMotifsGenome.pl ${PEAKS} ${ATYPE} ${OUTP}_motifs/ -size 50 -mask

# Deactivate environment
source deactivate
