# ATAC-Seq Tutorial

In this tutorial we will analyze the ATAC-Seq data from this [publication](https://www.ncbi.nlm.nih.gov/pubmed/27526324). The manuscript investigates the chromatin accessibility landscape of 13 human primary blood cell types. We will try to identify regulatory regions that govern hematopoietic differentiation.

## Data Set

The authors kindly provided a pre-processed [ATAC-Seq Count Matrix](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74912/suppl/GSE74912_ATACseq_All_Counts.txt.gz) which we will use in this tutorial. If you are interested in the basic Bioinformatics methods how ATAC-Seq data is aligned to the human reference genome, quality controlled and processed for peak calling and peak annotation you can have a look at our [ATAC-Seq pipeline](https://github.com/tobiasrausch/ATACseq). In this tutorial we will start with the raw count matrix of read counts per sample (columns) and peak (rows).

## Prepare the input count matrix

The [tutorial repository](https://github.com/tobiasrausch/ATACseq/tree/master/doc) contains a shell script to download the count matrix and subset the samples to one donor.

## Issues/Questions

For any question or error you may notice in the tutorial please open an issue in the [GitHub tutorial repository](https://github.com/tobiasrausch/ATACseq/tree/master/doc). Thanks.
