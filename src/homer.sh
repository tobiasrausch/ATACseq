#!/bin/bash

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <peaks.tsv> <align.bam> <genome.fa> <output prefix>"
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
ALIGN=${2}
HG=${3}
OUTP=${4}

# create tag directory (should we use normGC? only unique, keepOne?)
makeTagDirectory ${OUTP}/tagdir -genome ${HG} -checkGC ${ALIGN} 2> ${OUTP}.homer.log

# Annotated and normalized peaks
annotatePeaks.pl ${PEAKS} hg19 -size given -annStats ${OUTP}.homer.annStats -d ${OUTP}/tagdir > ${OUTP}.annotated.normalized 2>> ${OUTP}.homer.log

# create UCSC files
makeUCSCfile ${OUTP}/tagdir -style dnase -fsize 5e7 -o ${OUTP}.bedGraph 2>> ${OUTP}.homer.log.gz

# create IGV files
igvtools totdf ${OUTP}.bedGraph.gz ${OUTP}.tdf hg19

# Deactivate environment
source deactivate
