#!/bin/bash

if [ $# -ne 6 ]
then
    echo "**********************************************************************"
    echo "ATAC-Seq analysis pipeline."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <genome.fa> <outprefix> <control.rep1.final_peaks.narrowPeak> <control.rep2.final_peaks.narrowPeak> <treatment.rep1.final_peaks.narrowPeak> <treatment.rep2.final_peaks.narrowPeak>"
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PATH=${BASEDIR}/homer/bin:/g/funcgen/bin:${PATH}

# Custom parameters
THREADS=4

# CMD parameters
HG=${1}
OUTP=${2}
CP1=${3}
CP2=${4}
TP1=${5}
TP2=${6}
CB1=`echo ${CP1} | sed 's/_peaks.narrowPeak/.bam/'`
CB2=`echo ${CP2} | sed 's/_peaks.narrowPeak/.bam/'`
TB1=`echo ${TP1} | sed 's/_peaks.narrowPeak/.bam/'`
TB2=`echo ${TP2} | sed 's/_peaks.narrowPeak/.bam/'`
CT1=`echo ${CP1} | sed 's/\/\(.*\).pseudorep[12]_peaks.narrowPeak/\/\1\/tagdir\//' | sed 's/\/\(.*\).final_peaks.narrowPeak/\/\1\/tagdir\//'`
CT2=`echo ${CP2} | sed 's/\/\(.*\).pseudorep[12]_peaks.narrowPeak/\/\1\/tagdir\//' | sed 's/\/\(.*\).final_peaks.narrowPeak/\/\1\/tagdir\//'`
TT1=`echo ${TP1} | sed 's/\/\(.*\).pseudorep[12]_peaks.narrowPeak/\/\1\/tagdir\//' | sed 's/\/\(.*\).final_peaks.narrowPeak/\/\1\/tagdir\//'`
TT2=`echo ${TP2} | sed 's/\/\(.*\).pseudorep[12]_peaks.narrowPeak/\/\1\/tagdir\//' | sed 's/\/\(.*\).final_peaks.narrowPeak/\/\1\/tagdir\//'`
REPLICATES=1
if [ ${CT1} == ${CT2} ]
then
    CTAGDIR=${CT1}
    REPLICATES=0
else
    CTAGDIR="${CT1} ${CT2}"
fi
if [ ${TT1} == ${TT2} ]
then
    TTAGDIR=${TT1}
    REPLICATES=0
else
    TTAGDIR="${TT1} ${TT2}"
fi

echo ""
echo "Input Files and directories"
echo ""
ls ${CP1} ${CP2} ${TP1} ${TP2} ${CB1} ${CB2} ${TB1} ${TB2} ${CT1} ${CT2} ${TT1} ${TT2}
if [ $? -ne 0 ]
then
    echo ""
    echo "One of the input files does not exist!!!"
    echo ""
    exit -1
fi

# Programs
RSCR=${BASEDIR}/../R
PY3=${BASEDIR}/python3/bin/

# Tmp directory
if [ -n "${TMPDIR}" ]
then
    export TMP=${TMPDIR}
fi

# Merge filtered BAMs
samtools merge -@ ${THREADS} ${OUTP}.bam ${CB1} ${CB2} ${TB1} ${TB2}
samtools index -@ ${THREADS} ${OUTP}.bam

# Run stats
alfred qc -b ${BASEDIR}/../bed/tss.bed -r ${HG} -o ${OUTP}.bamStats.tsv.gz ${OUTP}.bam
MICOL=`zgrep "^ME" ${OUTP}.bamStats.tsv.gz | head -n 1 | tr '\t' '\n'  | awk '{print NR"\t"$0;}' | grep "MedianInsertSize" | cut -f 1`
ISIZE=`zgrep "^ME" ${OUTP}.bamStats.tsv.gz | tail -n 1 | tr '\t' '\n'  | awk '{print NR"\t"$0;}' | grep -P "^${MICOL}\t" | cut -f 2`
echo "Insert size" ${ISIZE}

# call peaks
macs2 callpeak -g hs --nomodel --keep-dup all -p 0.01 --shift 0 --extsize ${ISIZE} -n ${OUTP} -t ${OUTP}.bam

# run IDR
unset PYTHONPATH
export PATH=${PY3}:${PATH}
IDRTHRES=0.1
idr --samples ${CP1} ${CP2} --peak-list ${OUTP}_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}.control.idr --soft-idr-threshold ${IDRTHRES} --plot --use-best-multisummit-IDR --log-output-file ${OUTP}.control.idr.log
echo "Samples " ${CP1} ${CP2} >> ${OUTP}.control.idr.log
idr --samples ${TP1} ${TP2} --peak-list ${OUTP}_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}.treatment.idr --soft-idr-threshold ${IDRTHRES} --plot --use-best-multisummit-IDR --log-output-file ${OUTP}.treatment.idr.log
echo "Samples " ${TP1} ${TP2} >> ${OUTP}.treatment.idr.log

# filter peaks based on IDR
IDRCUT=`echo "-l(${IDRTHRES})/l(10)" | bc -l`
cat ${OUTP}.control.idr | awk '$12>='"${IDRCUT}"'' | cut -f 1-10 > ${OUTP}.control.peaks
cat ${OUTP}.treatment.idr | awk '$12>='"${IDRCUT}"'' | cut -f 1-10 > ${OUTP}.treatment.peaks

# estimate noise as #reads outside IDR peaks
cat ${OUTP}.treatment.idr ${OUTP}.control.idr | cut -f 1-3 | sort -k1,1V -k2,2n | uniq | awk '{print $1"\t"$2"\t"$3"\tPeak"NR;}' > ${OUTP}.idrpeaks.bed
alfred qc -b ${OUTP}.idrpeaks.bed -r ${HG} -o ${OUTP}.idrpeaks.gz ${OUTP}.bam

# subset peaks to IDR-filtered peaks
bedtools intersect -a ${OUTP}_peaks.narrowPeak -b ${OUTP}.idrpeaks.bed | sort -k1,1V -k2,2n | uniq > ${OUTP}.peaks

# quantify peaks (or -len 0 or -mask)
annotatePeaks.pl ${OUTP}.peaks hg19 -size given -noadj -raw -noann -nogene -d ${CTAGDIR} ${TTAGDIR} > ${OUTP}.peaks.quant

# normalized tag counts
annotatePeaks.pl ${OUTP}.peaks hg19 -size given -noann -nogene -d ${CTAGDIR} ${TTAGDIR} > ${OUTP}.peaks.normalized

# Annotated and normalized peaks
annotatePeaks.pl ${OUTP}.peaks hg19 -size given -annStats ${OUTP}.homer.annStats -d ${CTAGDIR} ${TTAGDIR} > ${OUTP}.annotated.normalized

# get differential peaks
if [ ${REPLICATES} -eq 1 ]; then
    # Replicates
    getDifferentialPeaksReplicates.pl -p ${OUTP}.peaks -genome hg19 -DESeq2 -t ${TTAGDIR} -b ${CTAGDIR} > ${OUTP}.up.differentialpeaks
    getDifferentialPeaksReplicates.pl -p ${OUTP}.peaks -genome hg19 -DESeq2 -t ${CTAGDIR} -b ${TTAGDIR} > ${OUTP}.down.differentialpeaks
else
    # No replicates
    getDifferentialPeaks ${OUTP}.peaks ${TT1} ${CT1} > ${OUTP}.up.differentialpeaks
    getDifferentialPeaks ${OUTP}.peaks ${CT1} ${TT1} > ${OUTP}.down.differentialpeaks
fi

# Annotate differential peaks
for DIRID in up down
do
    cut -f 1-4 ${OUTP}.${DIRID}.differentialpeaks | grep -v "^#" | sort -k2,2V -k3,3n | awk '{print $2"\t"$3"\t"$4"\t"$1;}' > ${OUTP}.${DIRID}.peaks
    annotatePeaks.pl ${OUTP}.${DIRID}.peaks hg19 -size given -annStats ${OUTP}.${DIRID}.homer.annStats -go ${OUTP}/go${DIRID} -d ${CTAGDIR} ${TTAGDIR} > ${OUTP}.${DIRID}.annotated.normalized.peaks
    mkdir -p ${OUTP}/motifs${DIRID}
    findMotifsGenome.pl ${OUTP}.${DIRID}.peaks hg19 ${OUTP}/motifs${DIRID}/ -size 50 -mask
done

