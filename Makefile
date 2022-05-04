SHELL := /bin/bash

# Targets
TARGETS = .conda .mamba .envs .homer
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/conda && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.mamba: .conda
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y -n base -c conda-forge mamba && touch .mamba

.envs: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -n atac freebayes=1.2.0 samtools=1.7 openssl=1.0 igvtools=2.3.93 alfred=0.1.17 idr=2.0.3 cutadapt=1.16 fastqc=0.11.7 bcftools=1.7 bowtie2=2.3.4.1 bbmap=37.90 biobambam=2.0.87 bedtools=2.27.1 jalview=2.10.3 datamash=1.1.0 ncurses=5.9 && mamba create -y -n atac2 macs2=2.1.1.20160309 homer=4.9.1 samtools=1.7 openssl=1.0 && touch .envs

.homer: .conda .mamba .envs
	export PATH=${PBASE}/conda/bin:${PATH} && source activate atac2 && perl ${PBASE}/conda/envs/atac2/share/homer-4.9.1-6/configureHomer.pl -install hg19 && perl ${PBASE}/conda/envs/atac2/share/homer-4.9.1-6/configureHomer.pl -install hg38 && source deactivate && touch .homer

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/
