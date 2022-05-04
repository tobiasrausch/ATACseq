#!/bin/bash

if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: $0 <align.bam> <output prefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Activate environment
export PATH=${BASEDIR}/../conda/bin:${PATH}
source activate atac

# Custom parameters
THREADS=4

# CMD parameters
INP=${1}
OUTP=${2}

# Only keep proper paired-ends
samtools sort -@ ${THREADS} -o ${OUTP}.namesort.bam -n ${INP}
samtools fixmate -@ ${THREADS} -r ${OUTP}.namesort.bam ${OUTP}.fixmate.bam
rm ${OUTP}.namesort.bam

# Generate pseudo-replicates
LRANDOM=`samtools idxstats ${INP} | awk '{SUM+=$3+$4;} END {print (int(SUM/4)+1);}'`
samtools view -@ ${THREADS} ${OUTP}.fixmate.bam |  sed 'N;s/\n/@\t@/' | shuf | split -d -l ${LRANDOM} - ${OUTP}.rep
rm ${OUTP}.fixmate.bam
samtools view -@ ${THREADS} -H ${INP} > ${OUTP}.rep1.sam
cat ${OUTP}.rep00 | sed 's/@\t@/\n/' >> ${OUTP}.rep1.sam
samtools view -@ ${THREADS} -b ${OUTP}.rep1.sam > ${OUTP}.rep1.bam
rm ${OUTP}.rep1.sam
samtools sort -@ ${THREADS} -o ${OUTP}.pseudorep1.bam ${OUTP}.rep1.bam
samtools index -@ ${THREADS} ${OUTP}.pseudorep1.bam
rm ${OUTP}.rep1.bam
samtools view -@ ${THREADS} -H ${INP} > ${OUTP}.rep2.sam
cat ${OUTP}.rep01 | sed 's/@\t@/\n/' >> ${OUTP}.rep2.sam
samtools view -@ ${THREADS} -b ${OUTP}.rep2.sam > ${OUTP}.rep2.bam
rm ${OUTP}.rep2.sam
samtools sort -@ ${THREADS} -o ${OUTP}.pseudorep2.bam ${OUTP}.rep2.bam
samtools index -@ ${THREADS} ${OUTP}.pseudorep2.bam
rm ${OUTP}.rep2.bam
rm ${OUTP}.rep00 ${OUTP}.rep01

# Deactivate environment
source deactivate
