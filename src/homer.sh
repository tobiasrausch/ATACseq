#!/bin/bash

if [ $# -ne 5 ]
then
    echo ""
    echo "Usage: $0 <hg19|mm10> <peaks.tsv> <align.bam> <genome.fa> <output prefix>"
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
ALIGN=${3}
HG=${4}
OUTP=${5}

# create tag directory (should we use normGC? only unique, keepOne?)
makeTagDirectory ${OUTP}/tagdir -genome ${HG} -checkGC ${ALIGN} 2> ${OUTP}.homer.log

# Annotated and normalized peaks
annotatePeaks.pl ${PEAKS} ${ATYPE} -size given -annStats ${OUTP}.homer.annStats -d ${OUTP}/tagdir > ${OUTP}.annotated.normalized 2>> ${OUTP}.homer.log

# Deactivate environment
source deactivate
