#!/bin/bash

if [ $# -lt 4 ]
then
    echo "Usage: $0 <genome.fa> <outprefix> -c <controlrep1.bam> -c <controlrep2.bam> ... -t <treatmentrep1.bam> -t <treatmentrep2.bam> ..."
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
if [ -n "${SCRATCHDIR}" ]
then
    export TMP=${SCRATCHDIR}
    echo "scratch directory" ${SCRATCHDIR}
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
CTD=""
LASTC=""
TREATMENT=""
TTD=""
LASTT=""
for ((  i = 0 ;  i < ${#OPTARR[@]};  i++  ))
do
    if [ ${OPTARR[$i]} == "-c" ]
    then
	TYPE="CONTROL"
    elif [ ${OPTARR[$i]} == "-t" ]
    then
	TYPE="TREATMENT"
    else
	TD=`echo ${OPTARR[$i]} | sed 's/\/[^\/]*$/\/tagdir\//'`
	if [ ${TYPE} == "CONTROL" ]
	then
	    CONTROL=${CONTROL}" "${OPTARR[$i]}
	    CTD=${CTD}" "${TD}
	    LASTC=${TD}
	else
	    TREATMENT=${TREATMENT}" "${OPTARR[$i]}
	    TTD=${TTD}" "${TD}
	    LASTT=${TD}
	fi
    fi
done

# Merge filtered BAMs
samtools merge ${OUTP}/${BAMID}.bam ${CONTROL} ${TREATMENT}
samtools index ${OUTP}/${BAMID}.bam

# Pseudo replicates for control
if [ `echo ${CONTROL} | sed 's/[ \t][ \t]*/\n/g' | wc -l | cut -f 1` -eq 1 ]
then
    LRANDOM=`samtools idxstats ${CONTROL} | awk '{SUM+=$3+$4;} END {print SUM/2}'`
    samtools view ${CONTROL} | shuf | split -d -l ${LRANDOM} - ${OUTP}/control.rep
    samtools view -H ${CONTROL} > ${OUTP}/control.rep1.sam
    cat ${OUTP}/control.rep00 >> ${OUTP}/control.rep1.sam
    samtools view -H ${CONTROL} > ${OUTP}/control.rep2.sam
    cat ${OUTP}/control.rep01 >> ${OUTP}/control.rep2.sam
    rm ${OUTP}/control.rep00 ${OUTP}/control.rep01
    samtools sort -o ${OUTP}/control.rep1.bam ${OUTP}/control.rep1.sam
    samtools index ${OUTP}/control.rep1.bam
    macs2 callpeak --gsize hs --nomodel --nolambda --keep-dup all --call-summits --name ${OUTP}/control.rep1 --treatment ${OUTP}/${BAMID}.bam
    samtools sort -o ${OUTP}/control.rep2.bam ${OUTP}/control.rep2.sam
    samtools index ${OUTP}/control.rep2.bam
    macs2 callpeak --gsize hs --nomodel --nolambda --keep-dup all --call-summits --name ${OUTP}/control.rep2 --treatment ${OUTP}/${BAMID}.bam
    rm ${OUTP}/control.rep1.sam ${OUTP}/control.rep2.sam
fi

# Pseudo replicates for treatment
if [ `echo ${TREATMENT} | sed 's/[ \t][ \t]*/\n/g' | wc -l | cut -f 1` -eq 1 ]
then
    LRANDOM=`samtools idxstats ${TREATMENT} | awk '{SUM+=$3+$4;} END {print SUM/2}'`
    samtools view ${TREATMENT} | shuf | split -d -l ${LRANDOM} - ${OUTP}/treatment.rep
    samtools view -H ${TREATMENT} > ${OUTP}/treatment.rep1.sam
    cat ${OUTP}/treatment.rep00 >> ${OUTP}/treatment.rep1.sam
    samtools view -H ${TREATMENT} > ${OUTP}/treatment.rep2.sam
    cat ${OUTP}/treatment.rep01 >> ${OUTP}/treatment.rep2.sam
    rm ${OUTP}/treatment.rep00 ${OUTP}/treatment.rep01
    samtools sort -o ${OUTP}/treatment.rep1.bam ${OUTP}/treatment.rep1.sam
    samtools index ${OUTP}/treatment.rep1.bam
    macs2 callpeak --gsize hs --nomodel --nolambda --keep-dup all --call-summits --name ${OUTP}/treatment.rep1 --treatment ${OUTP}/${BAMID}.bam
    samtools sort -o ${OUTP}/treatment.rep2.bam ${OUTP}/treatment.rep2.sam
    samtools index ${OUTP}/treatment.rep2.bam
    macs2 callpeak --gsize hs --nomodel --nolambda --keep-dup all --call-summits --name ${OUTP}/treatment.rep2 --treatment ${OUTP}/${BAMID}.bam
    rm ${OUTP}/treatment.rep1.sam ${OUTP}/treatment.rep2.sam
fi

# call peaks
macs2 callpeak --gsize hs --nomodel --nolambda --keep-dup all --call-summits --name ${OUTP}/${BAMID} --treatment ${OUTP}/${BAMID}.bam

# filter peaks
cd ${OUTP}
bedtools intersect -a ${BAMID}_peaks.narrowPeak -b <(zcat ${BASEDIR}/../bed/wgEncodeDacMapabilityConsensusExcludable.bed.gz) -wao | awk '$11=="."' | cut -f 1-10 | sort -k1,1V -k2,2n | uniq > ${BAMID}_peaks.narrowPeak.tmp && mv ${BAMID}_peaks.narrowPeak.tmp ${BAMID}_peaks.narrowPeak
cd ..

# filter peaks based on IDR
unset PYTHONPATH
export PATH=${PY3}:${PATH}
if [ `echo ${CONTROL} | sed 's/[ \t][ \t]*/\n/g' | wc -l | cut -f 1` -eq 1 ]
then
    idr --samples ${OUTP}/control.rep1_peaks.narrowPeak ${OUTP}/control.rep2_peaks.narrowPeak --peak-list ${OUTP}/${BAMID}_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}/${BAMID}.control.idr --plot --idr-threshold 0.01
else
    IDRSAMPLES=`echo ${CONTROL} | sed 's/.final.bam/_peaks.narrowPeak/g'`
    idr --samples ${IDRSAMPLES} --peak-list ${OUTP}/${BAMID}_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}/${BAMID}.control.idr --plot --idr-threshold 0.05
fi
if [ `echo ${TREATMENT} | sed 's/[ \t][ \t]*/\n/g' | wc -l | cut -f 1` -eq 1 ]
then
    idr --samples ${OUTP}/treatment.rep1_peaks.narrowPeak ${OUTP}/treatment.rep2_peaks.narrowPeak --peak-list ${OUTP}/${BAMID}_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}/${BAMID}.treatment.idr --plot --idr-threshold 0.01
else
    IDRSAMPLES=`echo ${TREATMENT} | sed 's/.final.bam/_peaks.narrowPeak/g'`
    idr --samples ${IDRSAMPLES} --peak-list ${OUTP}/${BAMID}_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}/${BAMID}.treatment.idr --plot --idr-threshold 0.05
fi
cat ${OUTP}/${BAMID}.treatment.idr ${OUTP}/${BAMID}.control.idr | cut -f 1-3 | sort -k1,1V -k2,2n | uniq > ${OUTP}/${BAMID}.idr.peaks
bedtools intersect -a ${OUTP}/${BAMID}_peaks.narrowPeak -b ${OUTP}/${BAMID}.idr.peaks -wao | awk '$11!="."' | cut -f 1-10 | sort -k1,1V -k2,2n | uniq > ${OUTP}/${BAMID}_peaks.narrowPeak.tmp && mv ${OUTP}/${BAMID}_peaks.narrowPeak.tmp ${OUTP}/${BAMID}_peaks.narrowPeak

# quantify peaks (or -len 0 or -mask)
annotatePeaks.pl ${OUTP}/${BAMID}_peaks.narrowPeak hg19 -size given -noadj -raw -noann -nogene -d ${CTD} ${TTD} > ${OUTP}/${BAMID}_peaks.narrowPeak.quant

# get differential peaks
getDifferentialPeaks ${OUTP}/${BAMID}_peaks.narrowPeak ${LASTT} ${LASTC} > ${OUTP}/${BAMID}.up.differentialpeaks
getDifferentialPeaks ${OUTP}/${BAMID}_peaks.narrowPeak ${LASTC} ${LASTT} > ${OUTP}/${BAMID}.down.differentialpeaks
getDifferentialPeaksReplicates.pl -p ${OUTP}/${BAMID}_peaks.narrowPeak -genome hg19 -DESeq2 -b ${CTD} -t ${TTD} > ${OUTP}/${BAMID}.up.deseq2
getDifferentialPeaksReplicates.pl -p ${OUTP}/${BAMID}_peaks.narrowPeak -genome hg19 -DESeq2 -b ${TTD} -t ${CTD} > ${OUTP}/${BAMID}.down.deseq2

# TF motif prediction
cd ${OUTP}
mkdir -p motifs
findMotifsGenome.pl ${BAMID}_peaks.narrowPeak hg19 motifs/ -size 50 -mask

# Clean-up tmp
if [ -n "${SCRATCHDIR}" ]
then
    ls ${SCRATCHDIR}
else
    rm -rf ${TMP}
fi
