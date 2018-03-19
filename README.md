ATAC-Seq Pipeline Installation
------------------------------

`git clone https://github.com/tobiasrausch/ATACseq.git`

`cd ATACseq`

`make all`


Building Promoter Regions for QC
--------------------------------

`cd bed/ && Rscript promoter.R`


Running the ATAC-Seq analysis pipeline for a single sample
----------------------------------------------------------

`./src/atac.sh <hg19|mm10> <read1.fq.gz> <read2.fq.gz> <genome.fa> <output prefix>`


Differential peak calling
-------------------------

Merge peaks across samples and create a raw count matrix.

`ls ./Sample1/Sample1.peaks ./Sample2/Sample2.peaks ./SampleN/SampleN.peaks > peaks.lst`

`ls ./Sample1/Sample1.bam ./Sample2/Sample2.bam ./SampleN/SampleN.bam > bams.lst`

`./src/count.sh peaks.lst bams.lst <output prefix>`

To call differential peaks on a count matrix, called counts.gz, using DESeq2 we first need to create a file with sample level information (sample.info). For instance, if you have 2 replicates per condition:

`echo -e "name\tcondition" > sample.info`

`zcat counts.gz | head -n 1 | cut -f 5- | tr '\t' '\n' | sed 's/.final$//' | awk '{print $0"\t"int((NR-1)/2);}' >> sample.info`

`Rscript R/dpeaks.R counts.gz sample.info`

