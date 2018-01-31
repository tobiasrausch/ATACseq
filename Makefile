SHELL := /bin/bash

# Targets
TARGETS = .homer .idr .igvtools
PBASE=$(shell pwd)

all:   	$(TARGETS)

.igvtools:
	cd src/ && wget 'http://data.broadinstitute.org/igv/projects/downloads/igvtools_2.3.91.zip' && unzip igvtools_2.3.91.zip && rm igvtools_2.3.91.zip && cd ../ && touch .igvtools

.homer:
	module load Perl && cd src/homer/ && perl configureHomer.pl -install homer && perl configureHomer.pl -install hg19 && cd ../../ && touch .homer

.idr:
	wget 'https://repo.continuum.io/archive/Anaconda3-4.3.0-Linux-x86_64.sh' && bash Anaconda3-4.3.0-Linux-x86_64.sh -b -p ${PBASE}/src/python3/ && rm Anaconda3-4.3.0-Linux-x86_64.sh && unset PYTHONPATH && source ${PBASE}/src/python3/bin/activate && cd src/ && wget 'https://github.com/nboley/idr/archive/2.0.3.zip' && unzip 2.0.3.zip && rm 2.0.3.zip && cd idr-2.0.3/ && python setup.py install && source deactivate && cd ../../ && touch .idr

clean:
	mv src/homer/configureHomer.pl . && rm -rf src/homer/ && mkdir -p src/homer/ && mv configureHomer.pl src/homer/
	rm -rf $(TARGETS) $(TARGETS:=.o) src/python3/ src/idr-2.0.3/ src/IGVTools/
