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
TARGETS = .fastqc .cutadapt .bowtie .picard .htslib .samtools .bcftools .bamStats

all:   	$(TARGETS)

.fastqc:
	cd src && wget 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip' && unzip fastqc_v0.11.5.zip && chmod 755 FastQC/fastqc && rm fastqc_v0.11.5.zip && cd ../ && touch .fastqc

.bamStats: $(BSTATSSOURCES)
	cd src/bamStats && make all && cd ../../ && touch .bamStats

.bowtie: $(BOWTIESOURCES)
	cd src/bowtie && make && cd ../../ && touch .bowtie

.cutadapt:
	mkdir src/cutadapt/ && virtualenv ${PBASE}/src/cutadapt/venv && ${PBASE}/src/cutadapt/venv/bin/pip install --install-option="--install-scripts=${PBASE}/src/cutadapt/exe" cutadapt==1.10 && touch .cutadapt

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
	cd src/modular-boost && ./b2 --clean-all
	rm -rf $(TARGETS) $(TARGETS:=.o) src/FastQC src/cutadapt/
