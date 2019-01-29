SHELL := /bin/bash

# Targets
TARGETS = .conda .channels .envs .homer
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/bin && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.channels: .conda
	export PATH=${PBASE}/bin/bin:${PATH} && conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda && touch .channels

.envs: .conda .channels
	export PATH=${PBASE}/bin/bin:${PATH} && conda create -y --prefix=${PBASE}/bin/envs/atac samtools=1.7 igvtools=2.3.93 alfred=0.1.17 idr=2.0.3 cutadapt=1.16 fastqc=0.11.7 bcftools=1.7 bowtie2=2.3.4.1 bbmap=37.90 biobambam=2.0.87 bedtools=2.27.1 jalview=2.10.3 datamash=1.1.0 && conda install -y --prefix=${PBASE}/bin/envs/atac -c conda-forge ncurses=5.9 && conda create -y --prefix=${PBASE}/bin/envs/atac2 macs2=2.1.1.20160309 homer=4.9.1 samtools=1.7 && touch .envs

.homer: .conda .channels .envs
	export PATH=${PBASE}/bin/bin:${PATH} && source activate ${PBASE}/bin/envs/atac2 && perl ${PBASE}/bin/envs/atac2/share/homer-4.9.1-6/configureHomer.pl -install hg19 && source deactivate && touch .homer

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) bin/
