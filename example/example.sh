#!/bin/bash

if [ $# -ne 1 ]
then
    echo "**********************************************************************"
    echo "Downloads one pair of FASTQ files from the below publication and runs"
    echo "it through the ATAC-Seq pipeline. The example does require a"
    echo "Bowtie2-indexed hg19 reference genome."
    echo ""
    echo "Publication: https://www.ncbi.nlm.nih.gov/pubmed/27526324"
    echo ""
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <genome.fa>"
    echo ""
    exit -1
fi

# Bowtie-index reference genome
HG=${1}

# Parse ENA study xml
wget -O PRJNA301969.xml 'http://www.ebi.ac.uk/ena/data/view/PRJNA301969&display=xml'
URL=`xmllint --xpath "//PROJECT/PROJECT_LINKS/PROJECT_LINK/XREF_LINK" PRJNA301969.xml | grep "ENA-FASTQ-FILES" -A 1 | grep "<ID>" | sed 's/^.*http/http/' | sed 's/fields=.*$/fields=fastq_ftp,sample_title/'`
rm PRJNA301969.xml

# Get fastq files
wget -O fastq.txt "${URL}"
cat fastq.txt | sed 's/ /_/' | sed 's/;/\t/' | sed 's/-/_/g' | tail -n +2 | head -n 1 > fastq.list
rm fastq.txt

# Download & run ATAC-Seq pipeline
if [ -f fastq.list ]
then
    if [ `cat fastq.list | wc -l | cut -f 1` -gt 0 ]
    then
	while read FQ1 FQ2 SAMPLE
	do
	    wget -O ${SAMPLE}.1.fq.gz "${FQ1}"
	    wget -O ${SAMPLE}.2.fq.gz "${FQ2}"
	    if [ -f ${SAMPLE}.1.fq.gz ]
	    then
		if [ -f ${SAMPLE}.2.fq.gz ]
		then
		    echo ${SAMPLE}
		    mkdir ${SAMPLE}
		    if [ -d ${SAMPLE} ]
		    then
			../src/atac.sh hg19 ${SAMPLE}.1.fq.gz ${SAMPLE}.2.fq.gz ${HG} ${SAMPLE}/${SAMPLE}
		    fi
		    rm ${SAMPLE}.2.fq.gz
		fi
		rm ${SAMPLE}.1.fq.gz
	    fi
	done < fastq.list
    fi
    rm fastq.list
fi
