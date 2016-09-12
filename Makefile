#VIRTUALENV=/g/software/bin/virtualenv
#PYTHON=/g/software/bin/python-2.7
#JAVA=/g/software/bin/java-8
VIRTUALENV=virtualenv
PYTHON=python
JAVA=java

# External sources
BOWTIESOURCES = $(wildcard src/bowtie/src/*.h)
PICARDSOURCES = $(wildcard src/picard/src/java/picard/*/*.java)
BSTATSSOURCES = $(wildcard src/bamStats/src/*.h)
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
SAMSOURCES = $(wildcard src/samtools/*.c) $(wildcard src/samtools/*.h)
BCFSOURCES = $(wildcard src/bcftools/*.c) $(wildcard src/bcftools/*.h)
BOOSTSOURCES = $(wildcard src/modular-boost/libs/iostreams/include/boost/iostreams/*.hpp)
PBASE=$(shell pwd)

# Targets
TARGETS = .fastqc .cutadapt .macs2 .bowtie .picard .htslib .samtools .bcftools .bamStats .java

all:   	$(TARGETS)

.java:
	mkdir src/java/ && cd src/java/ && ln -s ${JAVA} java && cd ../../ && touch .java

.fastqc:
	cd src && wget 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip' && unzip fastqc_v0.11.5.zip && chmod 755 FastQC/fastqc && rm fastqc_v0.11.5.zip && cd ../ && touch .fastqc

.bamStats: $(BSTATSSOURCES)
	cd src/bamStats && make all && cd ../../ && touch .bamStats

.bowtie: $(BOWTIESOURCES)
	cd src/bowtie && make && cd ../../ && touch .bowtie

.cutadapt:
	mkdir src/cutadapt/ && ${VIRTUALENV} -p ${PYTHON} ${PBASE}/src/venv && ${PBASE}/src/venv/bin/pip install --install-option="--install-scripts=${PBASE}/src/cutadapt" cutadapt==1.10 && touch .cutadapt

.macs2: .cutadapt
	. ${PBASE}/src/venv/bin/activate && pip install numpy && pip install Cython && pip install MACS2 && touch .macs2

.picard: $(PICARDSOURCES)
	cd src/picard && ./gradlew shadowJar && cd ../../ && touch .picard

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && make lib-static && cd ../../ && touch .htslib

.samtools: $(SAMSOURCES)
	cd src/samtools && make && cd ../../ && touch .samtools

.bcftools: $(BCFSOURCES)
	cd src/bcftools && make && cd ../../ && touch .bcftools

clean:
	cd src/picard && ./gradlew clean
	cd src/bowtie && make clean
	cd src/htslib && make clean
	cd src/samtools && make clean
	cd src/bcftools && make clean
	cd src/bamStats && make clean
	rm -rf $(TARGETS) $(TARGETS:=.o) src/FastQC src/cutadapt/ src/venv/ src/java/
