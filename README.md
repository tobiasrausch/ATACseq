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
