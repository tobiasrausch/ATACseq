#!/bin/bash

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <hg19|mm10> <genome.fa> <input.bam> <output prefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Activate environment
export PATH=${BASEDIR}/../bin/bin:${PATH}
source activate ${BASEDIR}/../bin/envs/atac

# Custom parameters
THREADS=4

# CMD parameters
ATYPE=${1}
HG=${2}
BAM=${3}
OUTP=${4}

# Call footprints
alfred tracks -c 2 -o ${OUTP}.footprint.bedGraph.gz ${BAM}
igvtools totdf ${OUTP}.footprint.bedGraph.gz ${OUTP}.footprint.tdf hg19

# call peaks on bedGraph
source deactivate
source activate ${BASEDIR}/../bin/envs/atac2

# Find suitable cutoff for > 10,000 peaks
macs2 bdgpeakcall --cutoff-analysis -g 75 -l 50 -i <(zcat ${OUTP}.footprint.bedGraph.gz) -o ${OUTP}.footprints_narrowPeak
CUTOFF=`tail -n +2 ${OUTP}.footprints_narrowPeak | awk '$2>10000' | head -n 1 | cut -f 1 | awk '{print int($1);}'`
echo ${CUTOFF}
# At least 5 spanning pairs
if [ ${CUTOFF} -lt 5 ]
then
    CUTOFF=5
fi

# Call peaks
macs2 bdgpeakcall -c ${CUTOFF} -g 75 -l 50 -i <(zcat ${OUTP}.footprint.bedGraph.gz) -o ${OUTP}.footprints_narrowPeak

# Get rid of the track line
tail -n +2 ${OUTP}.footprints_narrowPeak > ${OUTP}.footprints

# Annotate motifs in footprints
if [[ ${ATYPE} = *"hg"* ]]
then
    alfred annotate -m ${BASEDIR}/../motif/jaspar.gz -r ${HG} -o ${OUTP}.footprints.motifs -u ${OUTP}.motifs.footprints ${OUTP}.footprints
fi

source deactivate
