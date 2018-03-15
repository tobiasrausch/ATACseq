#!/bin/bash

if [ $# -ne 3 ]
then
    echo ""
    echo "Usage: $0 <replicate1.bam> <replicate2.bam> <output prefix>"
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
REP1=${1}
REP2=${2}
OUTP=${3}

# Merge BAMs
samtools merge -@ ${THREADS} ${OUTP}.merge.bam ${REP1} ${REP2}
samtools index -@ ${THREADS} ${OUTP}.merge.bam

# Estimate insert size
ISIZE=`samtools view -f 2 -F 3852 ${OUTP}.merge.bam | awk '$9>length($10)' | cut -f 9 | head -n 1000000 | awk '{SUM+=$1} END {print int(SUM/NR);}'`

# call peaks on replicates and merged BAM
source deactivate
source activate ${BASEDIR}/../bin/envs/atac2
for PEAKBAM in ${OUTP}.merge.bam ${REP1} ${REP2}
do
    PEAKN=${PEAKBAM}.suf
    macs2 callpeak -g hs --nomodel --keep-dup all -p 0.01 --shift 0 --extsize ${ISIZE} -n ${PEAKN} -t ${PEAKBAM}
    rm ${PEAKN}_summits.bed ${PEAKN}_peaks.xls
done
rm ${OUTP}.merge.bam ${OUTP}.merge.bam.bai
source deactivate
source activate ${BASEDIR}/../bin/envs/atac

# filter peaks based on IDR
IDRTHRES=0.1
idr --samples ${REP1}.suf_peaks.narrowPeak ${REP2}.suf_peaks.narrowPeak --peak-list ${OUTP}.merge.bam.suf_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}.idr --soft-idr-threshold ${IDRTHRES} --plot --use-best-multisummit-IDR --log-output-file ${OUTP}.idr.log
IDRCUT=`echo "-l(${IDRTHRES})/l(10)" | bc -l`
cat ${OUTP}.merge.bam.suf_peaks.narrowPeak | grep -w -Ff <(cat ${OUTP}.idr | awk '$12>='"${IDRCUT}"'' | cut -f 1-3) > ${OUTP}.peaks 
rm ${REP1}.suf_peaks.narrowPeak ${REP2}.suf_peaks.narrowPeak ${OUTP}.merge.bam.suf_peaks.narrowPeak
gzip ${OUTP}.idr

# Create UCSC track
echo "track type=narrowPeak visibility=3 db=hg19 name=\"${OUTP}\" description=\"${OUTP} narrowPeaks\"" | gzip -c > ${OUTP}.narrowPeak.ucsc.bed.gz
echo "browser position chr12:125400362-125403757" | gzip -c >> ${OUTP}.narrowPeak.ucsc.bed.gz
cat ${OUTP}.peaks | gzip -c >> ${OUTP}.narrowPeak.ucsc.bed.gz

# Deactivate environment
source deactivate
