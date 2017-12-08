#!/bin/bash

if [ $# -lt 4 ]
then
    echo "Usage: $0 <genome.fa> <outprefix> -c <controlrep1_folder> -c <controlrep2_folder> ... -t <treatmentrep1_folder> -t <treatmentrep2_folder> ..."
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PATH=${BASEDIR}/homer/bin:${PATH}

# Custom parameters
THREADS=4

# CMD parameters
HG=${1}
OUTP=${2}
shift
shift

# Programs
RSCR=${BASEDIR}/../R
PY3=${BASEDIR}/python3/bin/

# Tmp directory
DSTR=$(date +'%a_%y%m%d_%H%M')
if [ -n "${TMPDIR}" ]
then
    export TMP=${TMPDIR}
    echo "tmp directory" ${TMPDIR}
else
    export TMP=/tmp/tmp_atac_${DSTR}
    mkdir -p ${TMP}
fi

# Generate IDs
BAMID=`echo ${OUTP} | sed 's/$/.union/'`
mkdir -p ${OUTP}

# Get treatment and control BAMs
OPTARR=( $@ )
TYPE="CONTROL"
CONTROL=""
TREATMENT=""
for ((  i = 0 ;  i < ${#OPTARR[@]};  i++  ))
do
    if [ ${OPTARR[$i]} == "-c" ] ; then
	TYPE="CONTROL"
    elif [ ${OPTARR[$i]} == "-t" ] ; then
	TYPE="TREATMENT"
    else
	if [ ${TYPE} == "CONTROL" ] ; then
	    CONTROL=${CONTROL}" "`echo ${OPTARR[$i]} | sed 's/$/#/'`
	else
	    TREATMENT=${TREATMENT}" "`echo ${OPTARR[$i]} | sed 's/$/#/'`
	fi
    fi
done

# Merge filtered BAMs
samtools merge ${OUTP}/${BAMID}.bam `echo ${CONTROL} | sed 's/#/\/*.final.bam /g'` `echo ${TREATMENT} | sed 's/#/\/*.final.bam /g'`
samtools index ${OUTP}/${BAMID}.bam

# Run stats
alfred qc -b ${BASEDIR}/../bed/tss.bed -r ${HG} -o ${OUTP}/${OUTP}.bamStats ${OUTP}/${BAMID}.bam
Rscript ${RSCR}/isize.R ${OUTP}/${OUTP}.bamStats.isize.tsv
MICOL=`cat ${OUTP}/${OUTP}.bamStats.metrics.tsv | head -n 1 | tr '\t' '\n'  | awk '{print NR"\t"$0;}' | grep "MedianInsertSize" | cut -f 1`
ISIZE=`cat ${OUTP}/${OUTP}.bamStats.metrics.tsv | tail -n 1 | tr '\t' '\n'  | awk '{print NR"\t"$0;}' | grep -P "^${MICOL}\t" | cut -f 2`
rm ${OUTP}/${OUTP}.bamStats.coverage.tsv ${OUTP}/${OUTP}.bamStats.bedcov.tsv ${OUTP}/${OUTP}.bamStats.ontarget.tsv
rm ${OUTP}/${OUTP}.bamStats.mapq.tsv

# call peaks
macs2 callpeak -g hs --nomodel --keep-dup all -p 0.01 --shift 0 --extsize ${ISIZE} -n ${OUTP}/${BAMID} -t ${OUTP}/${BAMID}.bam

# run IDR
unset PYTHONPATH
export PATH=${PY3}:${PATH}
IDRTHRES=0.1
if [ `echo ${CONTROL} | tr '#' '\n' | grep "." | wc -l | cut -f 1` -eq 1 ] ; then
    IDRCONTROLSAMPLES=`echo ${CONTROL} | sed 's/#/\/*.pseudorep1_peaks.narrowPeak /'`" "`echo ${CONTROL} | sed 's/#/\/*.pseudorep2_peaks.narrowPeak /'`
else
    IDRCONTROLSAMPLES=`echo ${CONTROL} | sed 's/#/\/*.final_peaks.narrowPeak /g'`
fi
if [ `echo ${TREATMENT} | tr '#' '\n' | grep "." | wc -l | cut -f 1` -eq 1 ] ; then
    IDRTREATSAMPLES=`echo ${TREATMENT} | sed 's/#/\/*.pseudorep1_peaks.narrowPeak /'`" "`echo ${TREATMENT} | sed 's/#/\/*.pseudorep2_peaks.narrowPeak /'`
else
    IDRTREATSAMPLES=`echo ${TREATMENT} | sed 's/#/\/*.final_peaks.narrowPeak /g'`
fi
idr --samples ${IDRCONTROLSAMPLES} --peak-list ${OUTP}/${BAMID}_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}/${BAMID}.control.idr --soft-idr-threshold ${IDRTHRES} --plot --use-best-multisummit-IDR --log-output-file ${OUTP}/${BAMID}.control.idr.log
echo "Samples " ${IDRCONTROLSAMPLES} >> ${OUTP}/${BAMID}.control.idr.log
idr --samples ${IDRTREATSAMPLES} --peak-list ${OUTP}/${BAMID}_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}/${BAMID}.treatment.idr --soft-idr-threshold ${IDRTHRES} --plot --use-best-multisummit-IDR --log-output-file ${OUTP}/${BAMID}.treatment.idr.log
echo "Samples " ${IDRTREATSAMPLES} >> ${OUTP}/${BAMID}.treatment.idr.log

# filter peaks based on IDR
IDRCUT=`echo "-l(${IDRTHRES})/l(10)" | bc -l`
cat ${OUTP}/${BAMID}.control.idr | awk '$12>='"${IDRCUT}"'' | cut -f 1-10 > ${OUTP}/${BAMID}.control.peaks
cat ${OUTP}/${BAMID}.treatment.idr | awk '$12>='"${IDRCUT}"'' | cut -f 1-10 > ${OUTP}/${BAMID}.treatment.peaks

# estimate noise as #reads outside IDR peaks
cat ${OUTP}/${BAMID}.treatment.idr ${OUTP}/${BAMID}.control.idr | cut -f 1-3 | sort -k1,1V -k2,2n | uniq | awk '{print $1"\t"$2"\t"$3"\tPeak"NR;}' > ${OUTP}/${BAMID}.idrpeaks.bed
alfred qc -b ${OUTP}/${BAMID}.idrpeaks.bed -r ${HG} -o ${OUTP}/${OUTP}.idrpeaks ${OUTP}/${BAMID}.bam
rm ${OUTP}/${OUTP}.idrpeaks.coverage.tsv ${OUTP}/${OUTP}.idrpeaks.bedcov.tsv ${OUTP}/${OUTP}.idrpeaks.isize.tsv
rm ${OUTP}/${OUTP}.idrpeaks.mapq.tsv ${OUTP}/${OUTP}.idrpeaks.ontarget.tsv ${OUTP}/${OUTP}.idrpeaks.readlength.tsv

# subset peaks to IDR-filtered peaks
bedtools intersect -a ${OUTP}/${BAMID}_peaks.narrowPeak -b ${OUTP}/${BAMID}.idrpeaks.bed | sort -k1,1V -k2,2n | uniq > ${OUTP}/${BAMID}.peaks

# quantify peaks (or -len 0 or -mask)
annotatePeaks.pl ${OUTP}/${BAMID}.peaks hg19 -size given -noadj -raw -noann -nogene -d `echo ${CONTROL} | sed 's/#/\/tagdir\/ /g'` `echo ${TREATMENT} | sed 's/#/\/tagdir\/ /g'` > ${OUTP}/${BAMID}.peaks.quant

# normalized tag counts
annotatePeaks.pl ${OUTP}/${BAMID}.peaks hg19 -size given -noann -nogene -d `echo ${CONTROL} | sed 's/#/\/tagdir\/ /g'` `echo ${TREATMENT} | sed 's/#/\/tagdir\/ /g'` > ${OUTP}/${BAMID}.peaks.normalized

# Annotated and normalized peaks
annotatePeaks.pl ${OUTP}/${BAMID}.peaks hg19 -size given -annStats ${OUTP}/${BAMID}.${DIRID}.homer.annStats -d `echo ${CONTROL} | sed 's/#/\/tagdir\/ /g'` `echo ${TREATMENT} | sed 's/#/\/tagdir\/ /g'` > ${OUTP}/${BAMID}.annotated.normalized

# get differential peaks
if [ \( `echo ${CONTROL} | tr '#' '\n' | grep "." | wc -l | cut -f 1` -eq 1 \) -o \( `echo ${TREATMENT} | tr '#' '\n' | grep "." | wc -l | cut -f 1` -eq 1 \) ] ; then
    # No replicates
    getDifferentialPeaks ${OUTP}/${BAMID}.peaks `echo ${TREATMENT} | sed 's/#.*$/\/tagdir\/ /'` `echo ${CONTROL} | sed 's/#.*$/\/tagdir\/ /'` > ${OUTP}/${BAMID}.up.differentialpeaks
    getDifferentialPeaks ${OUTP}/${BAMID}.peaks `echo ${CONTROL} | sed 's/#.*$/\/tagdir\/ /'` `echo ${TREATMENT} | sed 's/#.*$/\/tagdir\/ /'` > ${OUTP}/${BAMID}.down.differentialpeaks
else
    # Replicates
    getDifferentialPeaksReplicates.pl -p ${OUTP}/${BAMID}.peaks -genome hg19 -DESeq2 -b `echo ${CONTROL} | sed 's/#/\/tagdir\/ /g'` -t `echo ${TREATMENT} | sed 's/#/\/tagdir\/ /g'` > ${OUTP}/${BAMID}.up.differentialpeaks
    getDifferentialPeaksReplicates.pl -p ${OUTP}/${BAMID}.peaks -genome hg19 -DESeq2 -b `echo ${TREATMENT} | sed 's/#/\/tagdir\/ /g'` -t `echo ${CONTROL} | sed 's/#/\/tagdir\/ /g'` > ${OUTP}/${BAMID}.down.differentialpeaks
fi

# Annotate differential peaks
for DIRID in up down
do
    cut -f 1-4 ${OUTP}/${BAMID}.${DIRID}.differentialpeaks | grep -v "^#" | sort -k2,2V -k3,3n | awk '{print $2"\t"$3"\t"$4"\t"$1;}' > ${OUTP}/${BAMID}.${DIRID}.peaks
    annotatePeaks.pl ${OUTP}/${BAMID}.${DIRID}.peaks hg19 -size given -annStats ${OUTP}/${BAMID}.${DIRID}.homer.annStats -go ${OUTP}/go${DIRID} -d `echo ${CONTROL} | sed 's/#/\/tagdir\/ /g'` `echo ${TREATMENT} | sed 's/#/\/tagdir\/ /g'` > ${OUTP}/${BAMID}.${DIRID}.annotated.normalized.peaks
    cd ${OUTP}
    mkdir -p motifs${DIRID}
    findMotifsGenome.pl ${BAMID}.${DIRID}.peaks hg19 motifs${DIRID}/ -size 50 -mask
    cd ../
done

# Clean-up tmp
if [ -n "${TMPDIR}" ]
then
    ls ${TMPDIR}
else
    rm -rf ${TMP}
fi
