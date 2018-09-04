#!/bin/bash

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <replicate1.bam> <replicate2.bam> <genome.fa> <output prefix>"
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
HG=${3}
OUTP=${4}

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
    #macs2 callpeak -g hs --nomodel --keep-dup all -p 0.01 --shift 0 --extsize ${ISIZE} -n ${PEAKN} -t ${PEAKBAM} -f BAMPE
    rm ${PEAKN}_summits.bed ${PEAKN}_peaks.xls
done
source deactivate
source activate ${BASEDIR}/../bin/envs/atac

# Saturated Peak Detection, significant peaks log2>=3 and -log10(p)>=3
PKTOTAL=`cat ${OUTP}.merge.bam.suf_peaks.narrowPeak | cut -f 1-3 | sort | uniq | wc -l | cut -f 1`
cat ${OUTP}.merge.bam.suf_peaks.narrowPeak | awk 'log($7)/log(2)>=2 && $8>=3' | cut -f 1-3 | sort -k1,1V -k2,2n | uniq > ${OUTP}.significant.peaks
cat ${REP1}.suf_peaks.narrowPeak | awk 'log($7)/log(2)>=1 && $8>=2' | cut -f 1-3 | sort -k1,1V -k2,2n | uniq > ${OUTP}.lenient.rep1.peaks
cat ${REP2}.suf_peaks.narrowPeak | awk 'log($7)/log(2)>=1 && $8>=2' | cut -f 1-3 | sort -k1,1V -k2,2n | uniq > ${OUTP}.lenient.rep2.peaks
RECALLREP1=`bedtools intersect -a ${OUTP}.significant.peaks -b ${OUTP}.lenient.rep1.peaks -wao | awk '$4!="."' | cut -f 1-3 | sort | uniq | wc -l | cut -f 1`
RECALLREP2=`bedtools intersect -a ${OUTP}.significant.peaks -b ${OUTP}.lenient.rep2.peaks -wao | awk '$4!="."' | cut -f 1-3 | sort | uniq | wc -l | cut -f 1`
SIGTOTAL=`cat ${OUTP}.significant.peaks | cut -f 1-3 | sort | uniq | wc -l | cut -f 1`
FRACREP1=`echo "${RECALLREP1} / ${SIGTOTAL}" | bc -l`
FRACREP2=`echo "${RECALLREP2} / ${SIGTOTAL}" | bc -l`
rm ${OUTP}.significant.peaks ${OUTP}.lenient.rep1.peaks ${OUTP}.lenient.rep2.peaks

# filter peaks based on IDR
IDRTHRES=0.1
idr --samples ${REP1}.suf_peaks.narrowPeak ${REP2}.suf_peaks.narrowPeak --peak-list ${OUTP}.merge.bam.suf_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}.idr --soft-idr-threshold ${IDRTHRES} --plot --use-best-multisummit-IDR --log-output-file ${OUTP}.idr.log
IDRCUT=`echo "-l(${IDRTHRES})/l(10)" | bc -l`
cat ${OUTP}.merge.bam.suf_peaks.narrowPeak | grep -w -Ff <(cat ${OUTP}.idr | awk '$12>='"${IDRCUT}"'' | cut -f 1-3) > ${OUTP}.peaks 
rm ${REP1}.suf_peaks.narrowPeak ${REP2}.suf_peaks.narrowPeak
mv ${OUTP}.merge.bam.suf_peaks.narrowPeak ${OUTP}.unfiltered.peaks
gzip ${OUTP}.unfiltered.peaks
gzip ${OUTP}.idr

# Fraction of reads in unfiltered peaks
alfred qc -b ${OUTP}.unfiltered.peaks.gz -r ${HG} -o ${OUTP}.bamStats.peaks.tsv.gz ${REP1}
FRACPEAK1=`zgrep "^ME" ${OUTP}.bamStats.peaks.tsv.gz  | datamash transpose | grep "^FractionInBed" | cut -f 2`
rm ${OUTP}.bamStats.peaks.tsv.gz
alfred qc -b ${OUTP}.unfiltered.peaks.gz -r ${HG} -o ${OUTP}.bamStats.peaks.tsv.gz ${REP2}
FRACPEAK2=`zgrep "^ME" ${OUTP}.bamStats.peaks.tsv.gz  | datamash transpose | grep "^FractionInBed" | cut -f 2`
rm ${OUTP}.merge.bam ${OUTP}.merge.bam.bai ${OUTP}.bamStats.peaks.tsv.gz

# Summarize peak statistics
echo -e "totpeaks\tfrip\tsigpeaks\trep1\trep2\trecallRep1\trecallRep2" > ${OUTP}.peaks.log
echo -e "${PKTOTAL}\t${FRACPEAK1},${FRACPEAK2}\t${SIGTOTAL}\t${RECALLREP1}\t${RECALLREP2}\t${FRACREP1}\t${FRACREP2}" >> ${OUTP}.peaks.log

# Create UCSC track
echo "track type=narrowPeak visibility=3 db=hg19 name=\"${OUTP}\" description=\"${OUTP} narrowPeaks\"" | gzip -c > ${OUTP}.narrowPeak.ucsc.bed.gz
echo "browser position chr12:125400362-125403757" | gzip -c >> ${OUTP}.narrowPeak.ucsc.bed.gz
cat ${OUTP}.peaks | gzip -c >> ${OUTP}.narrowPeak.ucsc.bed.gz

# Deactivate environment
source deactivate
