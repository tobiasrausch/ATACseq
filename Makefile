SHELL := /bin/bash

# Targets
TARGETS = .fastqc .homer .idr .picard 
PBASE=$(shell pwd)

all:   	$(TARGETS)

.fastqc:
	module load Java && cd src && wget 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip' && unzip fastqc_v0.11.5.zip && chmod 755 FastQC/fastqc && rm fastqc_v0.11.5.zip && cd ../ && touch .fastqc

.homer:
	module load Perl && cd src/homer/ && perl configureHomer.pl -install homer && perl configureHomer.pl -install hg19 && cd ../../ && touch .homer

.idr:
	wget 'https://repo.continuum.io/archive/Anaconda3-4.3.0-Linux-x86_64.sh' && bash Anaconda3-4.3.0-Linux-x86_64.sh -b -p ${PBASE}/src/python3/ && rm Anaconda3-4.3.0-Linux-x86_64.sh && unset PYTHONPATH && source ${PBASE}/src/python3/bin/activate && cd src/ && wget 'https://github.com/nboley/idr/archive/2.0.3.zip' && unzip 2.0.3.zip && rm 2.0.3.zip && cd idr-2.0.3/ && python setup.py install && source deactivate && cd ../../ && touch .idr

.picard:
	module load Java && mkdir -p src/picard/ && cd src/picard && wget -O picard.jar 'https://github.com/broadinstitute/picard/releases/download/2.8.3/picard.jar' && cd ../../ && touch .picard

clean:
	mv src/homer/configureHomer.pl . && rm -rf src/homer/ && mkdir -p src/homer/ && mv configureHomer.pl src/homer/
	rm -rf $(TARGETS) $(TARGETS:=.o) src/FastQC src/picard/ src/python3/ src/idr-2.0.3/
