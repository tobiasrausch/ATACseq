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

# call peaks
macs2 callpeak --gsize hs --nomodel --nolambda --keep-dup all --call-summits --name ${OUTP}/${BAMID} --treatment ${OUTP}/${BAMID}.bam
macs2 pileup --ifile ${OUTP}/${BAMID}.bam --ofile ${OUTP}/${BAMID}.bedGraph --format BAM --extsize 100

# extend peaks
bedtools slop -b 100 -i ${OUTP}/${BAMID}_summits.bed -g ${HG}.fai > ${OUTP}/${BAMID}.peaks

# filter peaks based on IDR
unset PYTHONPATH
export PATH=${PY3}:${PATH}
IDRSAMPLES=`echo ${CONTROL} | sed 's/.final.bam/_peaks.narrowPeak/g'`
idr --samples ${IDRSAMPLES} --peak-list ${OUTP}/${BAMID}_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}/${BAMID}.control.idr --plot --idr-threshold 0.05
IDRSAMPLES=`echo ${TREATMENT} | sed 's/.final.bam/_peaks.narrowPeak/g'`
idr --samples ${IDRSAMPLES} --peak-list ${OUTP}/${BAMID}_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}/${BAMID}.treatment.idr --plot --idr-threshold 0.05
cat ${OUTP}/${BAMID}.treatment.idr ${OUTP}/${BAMID}.control.idr | cut -f 1-3 | sort -k1,1V -k2,2n | uniq > ${OUTP}/${BAMID}.idr.peaks
bedtools intersect -a ${OUTP}/${BAMID}.peaks -b ${OUTP}/${BAMID}.idr.peaks | sort -k1,1V -k2,2n | uniq > ${OUTP}/${OUTP}.peaks.tmp && mv ${OUTP}/${OUTP}.peaks.tmp ${OUTP}/${BAMID}.peaks

# filter peaks against exclude regions
bedtools intersect -a ${OUTP}/${BAMID}.peaks -b <(zcat ${BASEDIR}/../bed/wgEncodeDacMapabilityConsensusExcludable.bed.gz) | cut -f 4 | sort | uniq > ${OUTP}/${OUTP}.remove
cat ${OUTP}/${BAMID}.peaks | grep -v -w -Ff ${OUTP}/${OUTP}.remove > ${OUTP}/${BAMID}.peaks.tmp && mv ${OUTP}/${BAMID}.peaks.tmp ${OUTP}/${BAMID}.peaks

# quantify peaks (or -len 0 or -mask)
annotatePeaks.pl ${OUTP}/${BAMID}.peaks hg19 -size given -noadj -raw -noann -nogene -d ${CTD} ${TTD} > ${OUTP}/${BAMID}.quant.raw.peaks

# get differential peaks
getDifferentialPeaks ${OUTP}/${BAMID}.peaks ${LASTT} ${LASTC} > ${OUTP}/${BAMID}.differentialpeaks
getDifferentialPeaksReplicates.pl -p ${OUTP}/${BAMID}.peaks -genome hg19 -DESeq2 -b ${CTD} -t ${TTD} > ${OUTP}/${BAMID}.deseq2

# TF motif prediction
cd ${OUTP}
mkdir -p motifs
findMotifsGenome.pl ${BAMID}.peaks hg19 motifs/ -size 50 -mask

# Clean-up tmp
if [ -n "${SCRATCHDIR}" ]
then
    ls ${SCRATCHDIR}
else
    rm -rf ${TMP}
fi
