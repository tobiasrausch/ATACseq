#!/bin/bash

if [ $# -ne 5 ]
then
    echo "**********************************************************************"
    echo "ATAC-Seq analysis pipeline."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Version: 0.1.1"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <hg19|mm10> <read1.fq.gz> <read2.fq.gz> <genome.fa> <output prefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# CMD parameters
ATYPE=${1}
FQ1=${2}
FQ2=${3}
HG=${4}
OUTP=${5}

# Align
${BASEDIR}/align.sh ${ATYPE} ${FQ1} ${FQ2} ${HG} ${OUTP}

# Generate pseudo-replicates
${BASEDIR}/pseudorep.sh ${OUTP}.final.bam ${OUTP}

# Call peaks and filter using IDR (replace pseudo-replicates with true biological replicates if available)
${BASEDIR}/peaks.sh ${OUTP}.pseudorep1.bam ${OUTP}.pseudorep2.bam ${HG} ${OUTP}

# Delete pseudo-replicates
rm ${OUTP}.pseudorep1.bam ${OUTP}.pseudorep1.bam.bai ${OUTP}.pseudorep2.bam ${OUTP}.pseudorep2.bam.bai

# Annotate peaks
${BASEDIR}/homer.sh ${OUTP}.peaks ${OUTP}.final.bam ${HG} ${OUTP}

# Motif discovery
${BASEDIR}/motif.sh ${ATYPE} ${OUTP}.peaks ${OUTP}

# Aggregate key QC metrics
${BASEDIR}/qc_globber.sh ${OUTP}
