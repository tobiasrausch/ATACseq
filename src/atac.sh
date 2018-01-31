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

# Add all required binaries
export PATH=${BASEDIR}/homer/bin:/g/funcgen/bin:${PATH}

# Custom parameters
THREADS=4
QUAL=30      # Mapping quality threshold

# CMD parameters
FQ1=${1}
FQ2=${2}
HG=${3}
OUTP=${4}

# Programs
RSCR=${BASEDIR}/../R
PY3=${BASEDIR}/python3/bin/
IGVTOOLS=${BASEDIR}/IGVTools/igvtools

# Tmp directory
if [ -n "${TMPDIR}" ]
then
    export TMP=${TMPDIR}
fi

# Generate IDs
FQ1ID=`echo ${OUTP} | sed 's/$/.fq1/'`
FQ2ID=`echo ${OUTP} | sed 's/$/.fq2/'`

# Adapter trimming
cutadapt -q 10 -m 15 -e 0.10 -a CTGTCTCTTATA -A CTGTCTCTTATA -o ${OUTP}.1.fq -p ${OUTP}.2.fq ${FQ1} ${FQ2} > ${OUTP}.cutadapt.log
gzip ${OUTP}.1.fq
gzip ${OUTP}.2.fq

# Bowtie
#bowtie2 --threads ${THREADS} --very-sensitive --maxins 2000  --no-discordant --no-mixed -x ${HG} -1 ${OUTP}.1.fq.gz -2 ${OUTP}.2.fq.gz
bowtie2 --threads ${THREADS} --local --maxins 2000 -x ${HG} -1 ${OUTP}.1.fq.gz -2 ${OUTP}.2.fq.gz | samtools view -bT ${HG} - > ${OUTP}.raw.bam

# Removed trimmed fastq
rm ${OUTP}.1.fq.gz ${OUTP}.2.fq.gz

# Sort & Index
samtools sort -@ ${THREADS} -o ${OUTP}.srt.bam ${OUTP}.raw.bam && rm ${OUTP}.raw.bam && samtools index ${OUTP}.srt.bam

# Mark duplicates
bammarkduplicates markthreads=${THREADS} tmpfile=${OUTP}_`date +'%H%M%S'` I=${OUTP}.srt.bam O=${OUTP}.rmdup.bam M=${OUTP}.metrics.tsv index=1 rmdup=0
rm ${OUTP}.srt.bam ${OUTP}.srt.bam.bai

# Run stats using unfiltered BAM
samtools idxstats ${OUTP}.rmdup.bam > ${OUTP}.idxstats
samtools flagstat ${OUTP}.rmdup.bam > ${OUTP}.flagstat

# Filter duplicates, mapping quality > QUAL, unmapped reads, chrM and unplaced contigs
CHRS=`cat ${BASEDIR}/../bed/tss.bed | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
samtools view -F 1804 -f 2 -q ${QUAL} -b ${OUTP}.rmdup.bam ${CHRS} > ${OUTP}.filt.bam
samtools index ${OUTP}.filt.bam
rm ${OUTP}.rmdup.bam ${OUTP}.rmdup.bam.bai

# Only keep proper paired-ends
samtools sort -o ${OUTP}.namesort.bam -n ${OUTP}.filt.bam
rm ${OUTP}.filt.bam ${OUTP}.filt.bam.bai
samtools fixmate -r ${OUTP}.namesort.bam ${OUTP}.fixmate.bam
rm ${OUTP}.namesort.bam
samtools sort -o ${OUTP}.final.bam ${OUTP}.fixmate.bam
samtools index ${OUTP}.final.bam

# Generate pseudo-replicates
LRANDOM=`samtools idxstats ${OUTP}.final.bam | awk '{SUM+=$3+$4;} END {print (int(SUM/4)+1);}'`
samtools view ${OUTP}.fixmate.bam |  sed 'N;s/\n/@\t@/' | shuf | split -d -l ${LRANDOM} - ${OUTP}.rep
rm ${OUTP}.fixmate.bam
samtools view -H ${OUTP}.final.bam > ${OUTP}.rep1.sam
cat ${OUTP}.rep00 | sed 's/@\t@/\n/' >> ${OUTP}.rep1.sam
samtools view -b ${OUTP}.rep1.sam > ${OUTP}.rep1.bam
rm ${OUTP}.rep1.sam
samtools sort -o ${OUTP}.pseudorep1.bam ${OUTP}.rep1.bam
samtools index ${OUTP}.pseudorep1.bam
rm ${OUTP}.rep1.bam
samtools view -H ${OUTP}.final.bam > ${OUTP}.rep2.sam
cat ${OUTP}.rep01 | sed 's/@\t@/\n/' >> ${OUTP}.rep2.sam
samtools view -b ${OUTP}.rep2.sam > ${OUTP}.rep2.bam
rm ${OUTP}.rep2.sam
samtools sort -o ${OUTP}.pseudorep2.bam ${OUTP}.rep2.bam
samtools index ${OUTP}.pseudorep2.bam
rm ${OUTP}.rep2.bam
rm ${OUTP}.rep00 ${OUTP}.rep01

# Run stats using filtered BAM, TSS enrichment, error rates, etc.
alfred qc -b ${BASEDIR}/../bed/tss.bed -r ${HG} -o ${OUTP}.bamStats.tsv.gz ${OUTP}.final.bam
MICOL=`zgrep "^ME" ${OUTP}.bamStats.tsv.gz | head -n 1 | tr '\t' '\n'  | awk '{print NR"\t"$0;}' | grep "MedianInsertSize" | cut -f 1`
ISIZE=`zgrep "^ME" ${OUTP}.bamStats.tsv.gz | tail -n 1 | tr '\t' '\n'  | awk '{print NR"\t"$0;}' | grep -P "^${MICOL}\t" | cut -f 2`

# FreeBayes
freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${HG} --genotype-qualities -b ${OUTP}.final.bam -v ${OUTP}.vcf
bgzip ${OUTP}.vcf
tabix ${OUTP}.vcf.gz

# Normalize VCF
bcftools norm -O z -o ${OUTP}.norm.vcf.gz -f ${GENOME} -m -both ${OUTP}.vcf.gz
tabix ${OUTP}.norm.vcf.gz
rm ${OUTP}.vcf.gz ${OUTP}.vcf.gz.tbi

# Fixed threshold filtering
bcftools filter -O z -o ${OUTP}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/AO<=2 || SAF<=2 || SAR<=2' ${OUTP}.norm.vcf.gz
tabix ${OUTP}.norm.filtered.vcf.gz
rm ${OUTP}.norm.vcf.gz ${OUTP}.norm.vcf.gz.tbi

# call peaks
for PEAKBAM in ${OUTP}.final.bam ${OUTP}.pseudorep1.bam ${OUTP}.pseudorep2.bam
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
idr --samples ${OUTP}.pseudorep1_peaks.narrowPeak ${OUTP}.pseudorep2_peaks.narrowPeak --peak-list ${OUTP}.final_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file ${OUTP}.idr --soft-idr-threshold ${IDRTHRES} --plot --use-best-multisummit-IDR --log-output-file ${OUTP}.idr.log
IDRCUT=`echo "-l(${IDRTHRES})/l(10)" | bc -l`
cat ${OUTP}.idr | awk '$12>='"${IDRCUT}"'' | cut -f 1-10 > ${OUTP}.peaks 

# estimate noise as #reads outside IDR peaks
cat ${OUTP}.peaks | awk '{print $1"\t"$2"\t"$3"\tPeak"NR;}' > ${OUTP}.idrpeaks.bed
alfred qc -b ${OUTP}.idrpeaks.bed -r ${HG} -o ${OUTP}.idrpeaks.gz ${OUTP}.final.bam

# annotate peaks using homer
cd ${OUTP}
annotatePeaks.pl ${OUTP}.peaks hg19 -annStats ${OUTP}.homer.annStats > ${OUTP}.annotated.peaks

# create tag directory (should we use normGC? only unique, keepOne?)
makeTagDirectory tagdir -genome ${HG} -checkGC ${OUTP}.final.bam
Rscript ${RSCR}/clonalTag.R tagdir/tagCountDistribution.txt
Rscript ${RSCR}/nuclfreq.R tagdir/tagFreq.txt
Rscript ${RSCR}/nuclfreq.R tagdir/tagFreqUniq.txt
Rscript ${RSCR}/autocor.R tagdir/tagAutocorrelation.txt
Rscript ${RSCR}/gc.R tagdir/genomeGCcontent.txt tagdir/tagGCcontent.txt

# Quantify peaks
annotatePeaks.pl ${OUTP}.peaks hg19 -size given -noann -nogene -d tagdir > ${OUTP}.peaks.normalized

# Annotated and normalized peaks
annotatePeaks.pl ${OUTP}.peaks hg19 -size given -annStats ${OUTP}.homer.annStats -d tagdir > ${OUTP}.annotated.normalized

# create UCSC files
makeUCSCfile tagdir -style dnase -fsize 5e7 -o ${OUTP}.bedGraph
echo "track type=narrowPeak visibility=3 db=hg19 name=\"${OUTP}\" description=\"${OUTP} narrowPeaks\"" | gzip -c > ${OUTP}.narrowPeak.ucsc.bed.gz
echo "browser position chr12:125400362-125403757" | gzip -c >> ${OUTP}.narrowPeak.ucsc.bed.gz
cat ${OUTP}.peaks | gzip -c >> ${OUTP}.narrowPeak.ucsc.bed.gz

# create IGV files
${IGVTOOLS} totdf ${OUTP}.bedGraph.gz ${OUTP}.tdf hg19

# TF motif prediction
mkdir -p motifs
findMotifsGenome.pl ${OUTP}.peaks hg19 motifs/ -size 50 -mask

# Clean-up tmp
if [ -n "${TMPDIR}" ]
then
    ls ${TMPDIR}
else
    rm -rf ${TMP}
fi
