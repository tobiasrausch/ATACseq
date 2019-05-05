ATAC-Seq Tutorial
-----------------

This is a tutorial to analyze an ATAC-Seq count matrix as produced by our [ATAC-Seq pipeline](https://github.com/tobiasrausch/ATACseq.git). The command requires [R Statistics](https://www.r-project.org/) and a number of additional R packages. To install R and the dependent packages just execute:

`make install`

If you already have R installed or want to install the packages manually just skip this step.

Downloading a published ATAC-Seq count matrix
---------------------------------------------

To download the tutorial data and subset to a small number of samples.

`./preprocess.sh`

This step produces the sample information table `blood.samples` and the actual count data `atac.data.gz`, which are required input files for the tutorial.


Building the ATAC-Seq Tutorial Book
-----------------------------------

You can build the ATAC-Seq book in HTML using

`make book`

You can then browse the tutorial content using

`firefox _book/index.html`

The online version of these HTML files you can find [here](https://tobiasrausch.com/)


Generating R scripts
--------------------

You can convert the Rmarkdown file to a standard Rscript file using

`make atac`


Credits
-------

This tutorial was created with contributions from [Markus Hsi-Yang Fritz](https://github.com/mhyfritz).
