#!/bin/bash

if [ $# -lt 3 ]
then
    echo ""
    echo "Usage: $0 <genome.fa> <output prefix> <sample1.bam> <sample2.bam> ... <sampleN.bam>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Activate environment
export PATH=${BASEDIR}/../conda/bin:${PATH}
source activate atac

# CMD parameters
GENOME=${1}
OUTP=${2}
shift 2

# Freebayes
freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities $@ -v ${OUTP}.vcf
cat ${OUTP}.vcf | grep "^#" > ${OUTP}.vcf.tmp
cat ${OUTP}.vcf| grep -v "^#" | sort -k1,1V -k2,2n >> ${OUTP}.vcf.tmp
mv ${OUTP}.vcf.tmp ${OUTP}.vcf
bgzip ${OUTP}.vcf
tabix ${OUTP}.vcf.gz

# Normalize VCF
bcftools norm -O z -o ${OUTP}.norm.vcf.gz -f ${GENOME} -m -both ${OUTP}.vcf.gz
tabix ${OUTP}.norm.vcf.gz
rm ${OUTP}.vcf.gz ${OUTP}.vcf.gz.tbi

# Filter VCF
bcftools filter -O z -o ${OUTP}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/INFO/AO<=2 || SAF<=2 || SAR<=2' ${OUTP}.norm.vcf.gz
tabix ${OUTP}.norm.filtered.vcf.gz
rm ${OUTP}.norm.vcf.gz ${OUTP}.norm.vcf.gz.tbi

# Deactivate environment
source deactivate
