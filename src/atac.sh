#!/bin/bash

if [ $# -ne 4 ]
then
    echo "**********************************************************************"
    echo "ATAC-Seq analysis pipeline."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <read1.fq.gz> <read2.fq.gz> <genome.fa> <output prefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PATH=${BASEDIR}/homer/bin:${PATH}

# Custom parameters
THREADS=4
QUAL=30      # Mapping quality threshold

# CMD parameters
FQ1=${1}
FQ2=${2}
HG=${3}
OUTP=${4}

# Programs
FREEBAYES=${BASEDIR}/freebayes/bin/freebayes
PICARD=${BASEDIR}/picard/picard.jar
FASTQC=${BASEDIR}/FastQC/fastqc
RSCR=${BASEDIR}/../R
PY3=${BASEDIR}/python3/bin/
IGVTOOLS=${BASEDIR}/IGVTools/igvtools

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
JAVAOPT="-Xms4g -Xmx8g -XX:ParallelGCThreads=${THREADS} -Djava.io.tmpdir=${TMP}"
PICARDOPT="MAX_RECORDS_IN_RAM=5000000 TMP_DIR=${TMP} VALIDATION_STRINGENCY=SILENT"

# Generate IDs
FQ1ID=`echo ${OUTP} | sed 's/$/.fq1/'`
FQ2ID=`echo ${OUTP} | sed 's/$/.fq2/'`
BAMID=`echo ${OUTP} | sed 's/$/.align/'`

# Fastqc
mkdir -p ${OUTP}/prefastqc/ && ${FASTQC} -t ${THREADS} -o ${OUTP}/prefastqc/ ${FQ1} && ${FASTQC} -t ${THREADS} -o ${OUTP}/prefastqc/ ${FQ2}

# Adapter trimming
cutadapt -q 10 -m 15 -e 0.10 -a CTGTCTCTTATA -A CTGTCTCTTATA -o ${OUTP}/${OUTP}.1.fq -p ${OUTP}/${OUTP}.2.fq ${FQ1} ${FQ2} > ${OUTP}/${OUTP}.cutadapt.log
gzip ${OUTP}/${OUTP}.1.fq
gzip ${OUTP}/${OUTP}.2.fq

# Fastqc
mkdir -p ${OUTP}/postfastqc/ && ${FASTQC} -t ${THREADS} -o ${OUTP}/postfastqc/ ${OUTP}/${OUTP}.1.fq.gz && ${FASTQC} -t ${THREADS} -o ${OUTP}/postfastqc/ ${OUTP}/${OUTP}.2.fq.gz

# Bowtie
#bowtie2 --threads ${THREADS} --very-sensitive --maxins 2000  --no-discordant --no-mixed -x ${HG} -1 ${OUTP}/${OUTP}.1.fq.gz -2 ${OUTP}/${OUTP}.2.fq.gz
bowtie2 --threads ${THREADS} --local --maxins 2000 -x ${HG} -1 ${OUTP}/${OUTP}.1.fq.gz -2 ${OUTP}/${OUTP}.2.fq.gz | samtools view -bT ${HG} - > ${OUTP}/${BAMID}.bam

# Removed trimmed fastq
rm ${OUTP}/${OUTP}.1.fq.gz ${OUTP}/${OUTP}.2.fq.gz

# Sort & Index
samtools sort -o ${OUTP}/${BAMID}.srt.bam ${OUTP}/${BAMID}.bam && rm ${OUTP}/${BAMID}.bam && samtools index ${OUTP}/${BAMID}.srt.bam

# Clean .bam file
java ${JAVAOPT} -jar ${PICARD} CleanSam I=${OUTP}/${BAMID}.srt.bam O=${OUTP}/${BAMID}.srt.clean.bam ${PICARDOPT} && rm ${OUTP}/${BAMID}.srt.bam*

# Mark duplicates
java ${JAVAOPT} -jar ${PICARD} MarkDuplicates I=${OUTP}/${BAMID}.srt.clean.bam O=${OUTP}/${BAMID}.srt.clean.rmdup.bam M=${OUTP}/${OUTP}.markdups.log PG=null MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ${PICARDOPT} && rm ${OUTP}/${BAMID}.srt.clean.bam* && samtools index ${OUTP}/${BAMID}.srt.clean.rmdup.bam

# Run stats using unfiltered BAM
samtools idxstats ${OUTP}/${BAMID}.srt.clean.rmdup.bam > ${OUTP}/${OUTP}.idxstats
samtools flagstat ${OUTP}/${BAMID}.srt.clean.rmdup.bam > ${OUTP}/${OUTP}.flagstat

# Filter duplicates, mapping quality > QUAL, unmapped reads, chrM and unplaced contigs
CHRS=`cat ${BASEDIR}/../bed/tss.bed | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
samtools view -F 1804 -f 2 -q ${QUAL} -b ${OUTP}/${BAMID}.srt.clean.rmdup.bam ${CHRS} > ${OUTP}/${BAMID}.filt.bam
samtools index ${OUTP}/${BAMID}.filt.bam
rm ${OUTP}/${BAMID}.srt.clean.rmdup.bam ${OUTP}/${BAMID}.srt.clean.rmdup.bam.bai

# Only keep proper paired-ends
samtools sort -o ${OUTP}/${BAMID}.namesort.bam -n ${OUTP}/${BAMID}.filt.bam
rm ${OUTP}/${BAMID}.filt.bam ${OUTP}/${BAMID}.filt.bam.bai
samtools fixmate -r ${OUTP}/${BAMID}.namesort.bam ${OUTP}/${BAMID}.fixmate.bam
rm ${OUTP}/${BAMID}.namesort.bam
samtools sort -o ${OUTP}/${BAMID}.final.bam ${OUTP}/${BAMID}.fixmate.bam
samtools index ${OUTP}/${BAMID}.final.bam

# Generate pseudo-replicates
LRANDOM=`samtools idxstats ${OUTP}/${BAMID}.final.bam | awk '{SUM+=$3+$4;} END {print (int(SUM/4)+1);}'`
samtools view ${OUTP}/${BAMID}.fixmate.bam |  sed 'N;s/\n/@\t@/' | shuf | split -d -l ${LRANDOM} - ${OUTP}/${BAMID}.rep
rm ${OUTP}/${BAMID}.fixmate.bam
samtools view -H ${OUTP}/${BAMID}.final.bam > ${OUTP}/${BAMID}.rep1.sam
cat ${OUTP}/${BAMID}.rep00 | sed 's/@\t@/\n/' >> ${OUTP}/${BAMID}.rep1.sam
samtools view -b ${OUTP}/${BAMID}.rep1.sam > ${OUTP}/${BAMID}.rep1.bam
rm ${OUTP}/${BAMID}.rep1.sam
samtools sort -o ${OUTP}/${BAMID}.pseudorep1.bam ${OUTP}/${BAMID}.rep1.bam
samtools index ${OUTP}/${BAMID}.pseudorep1.bam
rm ${OUTP}/${BAMID}.rep1.bam
samtools view -H ${OUTP}/${BAMID}.final.bam > ${OUTP}/${BAMID}.rep2.sam
cat ${OUTP}/${BAMID}.rep01 | sed 's/@\t@/\n/' >> ${OUTP}/${BAMID}.rep2.sam
samtools view -b ${OUTP}/${BAMID}.rep2.sam > ${OUTP}/${BAMID}.rep2.bam
rm ${OUTP}/${BAMID}.rep2.sam
samtools sort -o ${OUTP}/${BAMID}.pseudorep2.bam ${OUTP}/${BAMID}.rep2.bam
samtools index ${OUTP}/${BAMID}.pseudorep2.bam
rm ${OUTP}/${BAMID}.rep2.bam
rm ${OUTP}/${BAMID}.rep00 ${OUTP}/${BAMID}.rep01

# Run stats using filtered BAM, TSS enrichment, error rates, etc.
alfred -b ${BASEDIR}/../bed/tss.bed -r ${HG} -o ${OUTP}/${OUTP}.bamStats ${OUTP}/${BAMID}.final.bam
Rscript ${RSCR}/isize.R ${OUTP}/${OUTP}.bamStats.isize.tsv
MICOL=`cat ${OUTP}/${OUTP}.bamStats.metrics.tsv | head -n 1 | tr '\t' '\n'  | awk '{print NR"\t"$0;}' | grep "MedianInsertSize" | cut -f 1`
ISIZE=`cat ${OUTP}/${OUTP}.bamStats.metrics.tsv | tail -n 1 | tr '\t' '\n'  | awk '{print NR"\t"$0;}' | grep -P "^${MICOL}\t" | cut -f 2`

# FreeBayes
${FREEBAYES} --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${HG} --genotype-qualities -b ${OUTP}/${BAMID}.final.bam -v ${OUTP}/${BAMID}.vcf
bgzip ${OUTP}/${BAMID}.vcf
tabix ${OUTP}/${BAMID}.vcf.gz

# Normalize VCF
vt normalize ${OUTP}/${BAMID}.vcf.gz -r ${HG} | vt decompose_blocksub - | vt decompose - | vt uniq - | bgzip > ${OUTP}/${BAMID}.norm.vcf.gz
tabix ${OUTP}/${BAMID}.norm.vcf.gz
rm ${OUTP}/${BAMID}.vcf.gz ${OUTP}/${BAMID}.vcf.gz.tbi

# Fixed threshold filtering
bcftools filter -O z -o ${OUTP}/${BAMID}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/AO<=2 || SAF<=1 || SAR<=1' ${OUTP}/${BAMID}.norm.vcf.gz
tabix ${OUTP}/${BAMID}.norm.filtered.vcf.gz
rm ${OUTP}/${BAMID}.norm.vcf.gz ${OUTP}/${BAMID}.norm.vcf.gz.tbi

# call peaks
for PEAKBAM in ${OUTP}/${BAMID}.final.bam ${OUTP}/${BAMID}.pseudorep1.bam ${OUTP}/${BAMID}.pseudorep2.bam
do
    PEAKN=`echo ${PEAKBAM} | sed 's/.bam$//'`
    macs2 callpeak -g hs --nomodel --keep-dup all -p 0.01 --shift 0 --extsize ${ISIZE} -n ${PEAKN} -t ${PEAKBAM}

    # filter peaks
    bedtools intersect -v -a ${PEAKN}_peaks.narrowPeak -b <(zcat ${BASEDIR}/../bed/wgEncodeDacMapabilityConsensusExcludable.bed.gz) | sort -k1,1V -k2,2n | uniq > ${PEAKN}_peaks.narrowPeak.tmp && mv ${PEAKN}_peaks.narrowPeak.tmp ${PEAKN}_peaks.narrowPeak
done

# filter peaks based on IDR
unset PYTHONPATH
export PATH=${PY3}:${PATH}
IDRTHRES=0.1
idr --samples ${OUTP}/${BAMID}.pseudorep1_peaks.narrowPeak ${OUTP}/${BAMID}.pseudorep2_peaks.narrowPeak --peak-list ${OUTP}/${BAMID}.final_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}/${BAMID}.idr --soft-idr-threshold ${IDRTHRES} --plot --use-best-multisummit-IDR --log-output-file ${OUTP}/${BAMID}.idr.log
IDRCUT=`echo "-l(${IDRTHRES})/l(10)" | bc -l`
cat ${OUTP}/${BAMID}.idr | awk '$12>='"${IDRCUT}"'' | cut -f 1-10 > ${OUTP}/${BAMID}.peaks 

# annotate peaks using homer
cd ${OUTP}
annotatePeaks.pl ${BAMID}.peaks hg19 -annStats ${BAMID}.homer.annStats > ${BAMID}.annotated.peaks

# create tag directory (should we use normGC? only unique, keepOne?)
makeTagDirectory tagdir -genome ${HG} -checkGC ${BAMID}.final.bam
Rscript ${RSCR}/clonalTag.R tagdir/tagCountDistribution.txt
Rscript ${RSCR}/nuclfreq.R tagdir/tagFreq.txt
Rscript ${RSCR}/nuclfreq.R tagdir/tagFreqUniq.txt
Rscript ${RSCR}/autocor.R tagdir/tagAutocorrelation.txt
Rscript ${RSCR}/gc.R tagdir/genomeGCcontent.txt tagdir/tagGCcontent.txt

# create UCSC files
makeUCSCfile tagdir -style dnase -fsize 5e7 -o ${BAMID}.bedGraph
echo "track type=narrowPeak visibility=3 db=hg19 name=\"${BAMID}\" description=\"${BAMID} narrowPeaks\"" | gzip -c > ${BAMID}.narrowPeak.ucsc.bed.gz
echo "browser position chr12:125400362-125403757" | gzip -c >> ${BAMID}.narrowPeak.ucsc.bed.gz
cat ${BAMID}.peaks | gzip -c >> ${BAMID}.narrowPeak.ucsc.bed.gz

# create IGV files
${IGVTOOLS} totdf ${BAMID}.bedGraph.gz ${BAMID}.tdf hg19

# TF motif prediction
mkdir -p motifs
findMotifsGenome.pl ${BAMID}.peaks hg19 motifs/ -size 50 -mask

# Clean-up tmp
if [ -n "${SCRATCHDIR}" ]
then
    ls ${SCRATCHDIR}
else
    rm -rf ${TMP}
fi
