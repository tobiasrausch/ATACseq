#!/bin/bash

if [ $# -ne 1 ]
then
    echo ""
    echo "Usage: $0 <output prefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Parameters
OUTP=${1}

# Activate environment
export PATH=${BASEDIR}/../conda/bin:${PATH}
source activate atac

# Collect QC information
${BASEDIR}/qc_globber.py -p ${OUTP} > ${OUTP}.key.metrics

# Deactivate environment
source deactivate
