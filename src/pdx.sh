#!/bin/bash

if [ $# -lt 4 ]
then
    echo ""
    echo "Usage: $0 <human_mouse index folder> <outprefix> <read1.fq.gz> <read2.fq.gz>"
    echo ""
    exit -1
fi


SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PATH=/g/funcgen/bin/:${PATH}

# CMD parameters
HG=${1}
OP=${2}
FQ1=${3}
FQ2=${4}

# Index the human and mouse reference genomes (required only once)
# bbsplit.sh build=1 ref_x=hs37d5.fa ref_y=mm10.fa

# BBMap
module load BBMap
bbsplit.sh -Xmx64g usejni=t build=1 path=${HG} in=${FQ1} in2=${FQ2} basename=${OP}%_#.fq minratio=0.5 maxindel=100000 minhits=1 ambiguous2=all local

# Statistics
cat ${OP}x_1.fq | awk 'NR%4==1' > ${OP}x.reads
cat ${OP}y_1.fq | awk 'NR%4==1' > ${OP}y.reads
MOUSE=`sort ${OP}x.reads ${OP}x.reads ${OP}y.reads | uniq -u | wc -l | cut -f 1`
HUMAN=`sort ${OP}x.reads ${OP}y.reads ${OP}y.reads | uniq -u | wc -l | cut -f 1`
AMBIG=`sort ${OP}x.reads ${OP}y.reads | uniq -d | wc -l | cut -f 1`
echo "Unique mouse reads" ${MOUSE} > ${OP}.filter.stats
echo "Unique human reads" ${HUMAN} >> ${OP}.filter.stats
echo "Ambiguous reads" ${AMBIG} >> ${OP}.filter.stats
FRAC=`echo "${MOUSE} / ( ${MOUSE} + ${HUMAN} + ${AMBIG} )" | bc -l`
echo "Mouse fraction" ${FRAC} >> ${OP}.filter.stats
rm ${OP}x.reads ${OP}y.reads
rm ${OP}y_1.fq ${OP}y_2.fq

# gzip
mv ${OP}x_1.fq ${OP}.cleaned.1.fq
mv ${OP}x_2.fq ${OP}.cleaned.2.fq
gzip ${OP}.cleaned.1.fq ${OP}.cleaned.2.fq

