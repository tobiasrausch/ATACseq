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

Call differential peaks on that count matrix using DESeq2.
