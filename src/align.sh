#!/bin/bash

if [ $# -ne 5 ]
then
    echo ""
    echo "Usage: $0 <hg19|mm10> <read1.fq.gz> <read2.fq.gz> <genome.fa> <output prefix>"
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
QUAL=30      # Mapping quality threshold

# CMD parameters
ATYPE=${1}
FQ1=${2}
FQ2=${3}
HG=${4}
OUTP=${5}

# Generate IDs
FQ1ID=`echo ${OUTP} | sed 's/$/.fq1/'`
FQ2ID=`echo ${OUTP} | sed 's/$/.fq2/'`

# Fastqc
mkdir -p ${OUTP}_prefastqc/ && fastqc -t ${THREADS} -o ${OUTP}_prefastqc/ ${FQ1} && fastqc -t ${THREADS} -o ${OUTP}_prefastqc/ ${FQ2}

# Adapter trimming
cutadapt -q 10 -m 15 -e 0.10 -a CTGTCTCTTATA -A CTGTCTCTTATA -o ${OUTP}.1.fq -p ${OUTP}.2.fq ${FQ1} ${FQ2} | gzip -c > ${OUTP}.cutadapt.log.gz
gzip ${OUTP}.1.fq
gzip ${OUTP}.2.fq

# Fastqc
mkdir -p ${OUTP}_postfastqc/ && fastqc -t ${THREADS} -o ${OUTP}_postfastqc/ ${OUTP}.1.fq.gz && fastqc -t ${THREADS} -o ${OUTP}_postfastqc/ ${OUTP}.2.fq.gz

# Bowtie
#bowtie2 --threads ${THREADS} --very-sensitive --maxins 2000  --no-discordant --no-mixed -x ${HG} -1 ${OUTP}.1.fq.gz -2 ${OUTP}.2.fq.gz
bowtie2 --threads ${THREADS} --local --maxins 2000 -x ${HG} -1 ${OUTP}.1.fq.gz -2 ${OUTP}.2.fq.gz 2> ${OUTP}.bowtie.log | samtools view -bT ${HG} - > ${OUTP}.raw.bam

# Removed trimmed fastq
rm ${OUTP}.1.fq.gz ${OUTP}.2.fq.gz

# Sort & Index
samtools sort -@ ${THREADS} -o ${OUTP}.srt.bam ${OUTP}.raw.bam && rm ${OUTP}.raw.bam && samtools index -@ ${THREADS} ${OUTP}.srt.bam

# Mark duplicates
bammarkduplicates markthreads=${THREADS} tmpfile=${OUTP}_`date +'%H%M%S'` I=${OUTP}.srt.bam O=${OUTP}.rmdup.bam M=${OUTP}.rmdup.log index=1 rmdup=0
rm ${OUTP}.srt.bam ${OUTP}.srt.bam.bai

# Run stats using unfiltered BAM
samtools idxstats ${OUTP}.rmdup.bam > ${OUTP}.idxstats
samtools flagstat -@ ${THREADS} ${OUTP}.rmdup.bam > ${OUTP}.flagstat

# Run stats on unfiltered BAM
alfred qc -r ${HG} -o ${OUTP}.bamStats.unfiltered.tsv.gz ${OUTP}.rmdup.bam

# Filter duplicates, mapping quality > QUAL, unmapped reads, chrM and unplaced contigs
CHRS=`zcat ${BASEDIR}/../bed/${ATYPE}.promoter.bed.gz | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
samtools view -@ ${THREADS} -F 1804 -f 2 -q ${QUAL} -b ${OUTP}.rmdup.bam ${CHRS} > ${OUTP}.final.bam
samtools index -@ ${THREADS} ${OUTP}.final.bam
rm ${OUTP}.rmdup.bam ${OUTP}.rmdup.bam.bai

# Run stats using filtered BAM using promoter regions
alfred qc -b ${BASEDIR}/../bed/hg19.promoter.bed.gz -r ${HG} -o ${OUTP}.bamStats.promoters.tsv.gz ${OUTP}.final.bam

# Create browser tracks (ToDo)
alfred tracks -o ${OUTP}.bedGraph.gz ${OUTP}.final.bam
igvtools totdf ${OUTP}.bedGraph.gz ${OUTP}.tdf ${ATYPE}

# Deactivate environment
source deactivate


