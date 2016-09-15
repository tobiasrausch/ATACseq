#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Usage: $0 <read1.fq.gz> <read2.fq.gz> <genome.fa> <outprefix>"
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/g/software/linux/pack/python-2.7/lib/:/usr/lib64/
export PATH=${BASEDIR}/homer/bin:${BASEDIR}/gs/bin:${BASEDIR}/blat:${BASEDIR}/weblogo:${PATH}

# Custom parameters
THREADS=4
QUAL=10      # Mapping quality threshold

# CMD parameters
FQ1=${1}
FQ2=${2}
HG=${3}
OUTP=${4}

# Programs
PICARD=${BASEDIR}/picard/build/libs/picard.jar
SAM=${BASEDIR}/samtools/samtools
BCF=${BASEDIR}/bcftools/bcftools
FASTQC=${BASEDIR}/FastQC/fastqc
BOWTIE=${BASEDIR}/bowtie/bowtie2
BEDTOOLS=${BASEDIR}/bedtools/bin/bedtools
CUTADAPT=${BASEDIR}/cutadapt/cutadapt
BAMSTATS=${BASEDIR}/bamStats/src/bamStats
BAMSTATR=${BASEDIR}/bamStats/R
ACT=${BASEDIR}/venv/bin/activate
JAVA=${BASEDIR}/java/java

# Tmp directory
DSTR=$(date +'%a_%y%m%d_%H%M')
export TMP=/tmp/tmp_atac_${DSTR}
mkdir -p ${TMP}
JAVAOPT="-Xms4g -Xmx8g -XX:ParallelGCThreads=${THREADS} -Djava.io.tmpdir=${TMP}"
PICARDOPT="MAX_RECORDS_IN_RAM=5000000 TMP_DIR=${TMP} VALIDATION_STRINGENCY=SILENT"

# Generate IDs
FQ1ID=`echo ${OUTP} | sed 's/$/.fq1/'`
FQ2ID=`echo ${OUTP} | sed 's/$/.fq2/'`
BAMID=`echo ${OUTP} | sed 's/$/.align/'`

# Fastqc
mkdir -p ${OUTP}/prefastqc/ && ${FASTQC} -t ${THREADS} -o ${OUTP}/prefastqc/ ${FQ1} && ${FASTQC} -t ${THREADS} -o ${OUTP}/prefastqc/ ${FQ2}

# Adapter trimming
${CUTADAPT} --quiet -q 20 -m 35 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -o ${OUTP}/${OUTP}.1.fq -p ${OUTP}/${OUTP}.2.fq ${FQ1} ${FQ2} && gzip ${OUTP}/${OUTP}.1.fq && gzip ${OUTP}/${OUTP}.2.fq

# Fastqc
mkdir -p ${OUTP}/postfastqc/ && ${FASTQC} -t ${THREADS} -o ${OUTP}/postfastqc/ ${OUTP}/${OUTP}.1.fq.gz && ${FASTQC} -t ${THREADS} -o ${OUTP}/postfastqc/ ${OUTP}/${OUTP}.2.fq.gz

# Bowtie
${BOWTIE} --threads ${THREADS} --very-sensitive --maxins 2000  --no-discordant --no-mixed -x ${HG} -1 ${OUTP}/${OUTP}.1.fq.gz -2 ${OUTP}/${OUTP}.2.fq.gz | samtools view -bT ${HG} - > ${OUTP}/${BAMID}.bam

# Removed trimmed fastq
rm ${OUTP}/${OUTP}.1.fq.gz ${OUTP}/${OUTP}.2.fq.gz

# Sort & Index
${SAM} sort -o ${OUTP}/${BAMID}.srt.bam ${OUTP}/${BAMID}.bam && rm ${OUTP}/${BAMID}.bam && ${SAM} index ${OUTP}/${BAMID}.srt.bam

# Clean .bam file
${JAVA} ${JAVAOPT} -jar ${PICARD} CleanSam I=${OUTP}/${BAMID}.srt.bam O=${OUTP}/${BAMID}.srt.clean.bam ${PICARDOPT} && rm ${OUTP}/${BAMID}.srt.bam*

# Mark duplicates
${JAVA} ${JAVAOPT} -jar ${PICARD} MarkDuplicates I=${OUTP}/${BAMID}.srt.clean.bam O=${OUTP}/${BAMID}.srt.clean.rmdup.bam M=${OUTP}/${OUTP}.markdups.log PG=null MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ${PICARDOPT} && rm ${OUTP}/${BAMID}.srt.clean.bam* && ${SAM} index ${OUTP}/${BAMID}.srt.clean.rmdup.bam

# Run stats using unfiltered BAM
${SAM} idxstats ${OUTP}/${BAMID}.srt.clean.rmdup.bam > ${OUTP}/${OUTP}.idxstats
${SAM} flagstat ${OUTP}/${BAMID}.srt.clean.rmdup.bam > ${OUTP}/${OUTP}.flagstat

# Filter duplicates, mapping quality > QUAL, unmapped reads, chrM and unplaced contigs
CHRS=`cat ${BASEDIR}/../bed/tss.bed | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
${SAM} view -F 1024 -q ${QUAL} -b ${OUTP}/${BAMID}.srt.clean.rmdup.bam ${CHRS} > ${OUTP}/${BAMID}.final.bam
${SAM} index ${OUTP}/${BAMID}.final.bam
rm ${OUTP}/${BAMID}.srt.clean.rmdup.bam ${OUTP}/${BAMID}.srt.clean.rmdup.bam.bai

# Run stats using filtered BAM, TSS enrichment, error rates, etc.
${BAMSTATS} -b ${BASEDIR}/../bed/tss.bed -r ${HG} -o ${OUTP}/${OUTP}.bamStats ${OUTP}/${BAMID}.final.bam
Rscript ${BAMSTATR}/isize.R ${OUTP}/${OUTP}.bamStats.isize.tsv
Rscript ${BAMSTATR}/mapq.R ${OUTP}/${OUTP}.bamStats.mapq.tsv

# call peaks
source ${ACT}
#macs2 callpeak --gsize hs --nomodel --shift -100 --extsize 200 --broad --name ${OUTP}/${BAMID} --treatment ${OUTP}/${BAMID}.final.bam
macs2 callpeak --gsize hs --nomodel --nolambda --keep-dup all --call-summits --name ${OUTP}/${BAMID} --treatment ${OUTP}/${BAMID}.final.bam
macs2 pileup --ifile ${OUTP}/${BAMID}.final.bam --ofile ${OUTP}/${BAMID}.bedGraph --format BAM --extsize 100

# extend peaks
cd ${OUTP}
${BEDTOOLS} slop -b 100 -i ${BAMID}_summits.bed -g ${HG}.fai > ${BAMID}.peaks

# filter peaks
${BEDTOOLS} intersect -a ${BAMID}.peaks -b <(zcat ${BASEDIR}/../bed/wgEncodeDacMapabilityConsensusExcludable.bed.gz) | cut -f 4 | sort | uniq > ${OUTP}.remove
cat ${BAMID}.peaks | grep -v -w -Ff ${OUTP}.remove > ${BAMID}.peaks.tmp && mv ${BAMID}.peaks.tmp ${BAMID}.peaks

# annotate peaks using homer
annotatePeaks.pl ${BAMID}.peaks hg19 -annStats ${BAMID}.homer.annStats > ${BAMID}.annotated.peaks

# TF motif prediction
mkdir -p motifs
findMotifsGenome.pl ${BAMID}.peaks hg19 motifs/ -size 50 -mask

# Clean-up tmp
rm -rf ${TMP}
