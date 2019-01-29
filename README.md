ATAC-Seq Pipeline Installation
------------------------------

`git clone https://github.com/tobiasrausch/ATACseq.git`

`cd ATACseq`

`make all`

If one of the above commands fail your operating system probably lacks some build essentials. These are usually pre-installed but if you lack them you need to install these. For instance, for Ubuntu this would require:

`apt-get install build-essential g++ git wget unzip`


Building promoter regions for QC and downloading motifs
-------------------------------------------------------

To annotate motifs and estimate TSS enrichments some simple scripts are included in this repository to download these databases.

`cd bed/ && Rscript promoter.R && cd ..`

`cd motif/ && ./downloadMotifs.sh && cd ..`


Running the ATAC-Seq analysis pipeline for a single sample
----------------------------------------------------------

`./src/atac.sh <hg19|mm10> <read1.fq.gz> <read2.fq.gz> <genome.fa> <output prefix>`


Plotting the key ATAC-Seq Quality Control metrics
-------------------------------------------------

If you have multiple output folders (one for each ATAC-Seq sample) you can simply concatenate the QC metrics of each sample.

`head -n 1 ./*/*.key.metrics | grep "TssEnrichment" | uniq > summary.tsv`

`cat ./*/*.key.metrics | grep -v "TssEnrichment" >> summary.tsv`

To plot the distribution for all QC parameters.

`Rscript R/metrics.R summary.tsv`


Differential peak calling
-------------------------

Merge peaks across samples and create a raw count matrix.

`ls ./Sample1/Sample1.peaks ./Sample2/Sample2.peaks ./SampleN/SampleN.peaks > peaks.lst`

`ls ./Sample1/Sample1.bam ./Sample2/Sample2.bam ./SampleN/SampleN.bam > bams.lst`

`./src/count.sh hg19 peaks.lst bams.lst <output prefix>`

To call differential peaks on a count matrix for TSS peaks, called counts.tss.gz, using DESeq2 we first need to create a file with sample level information (sample.info). For instance, if you have 2 replicates per condition:

`echo -e "name\tcondition" > sample.info`

`zcat counts.tss.gz | head -n 1 | cut -f 5- | tr '\t' '\n' | sed 's/.final$//' | awk '{print $0"\t"int((NR-1)/2);}' >> sample.info`

`Rscript R/dpeaks.R counts.tss.gz sample.info`


Footprinting
------------

The pipeline also has a module to call footprints of nucleosome occupancy or transcription factor binding occupancy that are annotated for motif hits.

`./src/footprints.sh <hg19|mm10> <genome.fa> <input.bam> <output prefix>`


Citation
--------

Tobias Rausch, Markus Hsi-Yang Fritz, Jan O Korbel, Vladimir Benes.       
[Alfred: Interactive multi-sample BAM alignment statistics, feature counting and feature annotation for long- and short-read sequencing.](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty1007/5232224)      
Bioinformatics.


License
-------
This ATAC-Seq pipeline is distributed under the GPLv3.