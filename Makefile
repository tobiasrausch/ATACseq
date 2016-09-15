#VIRTUALENV=/g/software/bin/virtualenv
#PYTHON=/g/software/bin/python-2.7
#JAVA=/g/software/bin/java-8
VIRTUALENV=virtualenv
PYTHON=python
JAVA=java

# External sources
BOWTIESOURCES = $(wildcard src/bowtie/src/*.h)
BEDSOURCES = $(wildcard src/bedtools/src/*/*.cpp)
PICARDSOURCES = $(wildcard src/picard/src/java/picard/*/*.java)
BSTATSSOURCES = $(wildcard src/bamStats/src/*.h)
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
SAMSOURCES = $(wildcard src/samtools/*.c) $(wildcard src/samtools/*.h)
BCFSOURCES = $(wildcard src/bcftools/*.c) $(wildcard src/bcftools/*.h)
BOOSTSOURCES = $(wildcard src/modular-boost/libs/iostreams/include/boost/iostreams/*.hpp)
PBASE=$(shell pwd)

# Targets
TARGETS = .fastqc .cutadapt .macs2 .bedtools .gs .weblogo .blat .homer .bowtie .picard .htslib .samtools .bcftools .bamStats .java

all:   	$(TARGETS)

.java:
	mkdir src/java/ && cd src/java/ && ln -s ${JAVA} java && cd ../../ && touch .java

.fastqc:
	cd src && wget 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip' && unzip fastqc_v0.11.5.zip && chmod 755 FastQC/fastqc && rm fastqc_v0.11.5.zip && cd ../ && touch .fastqc

.bamStats: $(BSTATSSOURCES)
	cd src/bamStats && make all && cd ../../ && touch .bamStats

.bedtools: $(BEDSOURCES)
	cd src/bedtools && make all && cd ../../ && touch .bedtools

.gs:
	cd src/ && wget wget -O ghostscript-9.19.tar.gz 'https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs919/ghostscript-9.19.tar.gz' && tar -xzf ghostscript-9.19.tar.gz && rm ghostscript-9.19.tar.gz && cd ghostscript-9.19/ && ./configure --prefix=${PBASE}/src/gs/ && make && make install && cd ../ && rm -rf ghostscript-9.19/ && cd ../ && touch .gs

.weblogo: .gs
	cd src/ && wget 'http://weblogo.berkeley.edu/release/weblogo.2.8.2.tar.gz' && tar -xzf weblogo.2.8.2.tar.gz && rm -rf weblogo.2.8.2.tar.gz && cd ../ && touch .weblogo

.blat: .blat
	mkdir src/blat/ && cd src/blat/ && wget 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat' && chmod 755 blat && cd ../../ && touch .blat

.homer: .gs .weblogo .blat
	export PATH=${PBASE}/src/gs/bin:${PBASE}/src/weblogo:${PBASE}/src/blat:${PATH} && cd src/homer/ && perl configureHomer.pl -install homer && perl configureHomer.pl -install hg19 && cd ../../ && touch .homer

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
	cd src/bedtools && make clean
	cd src/htslib && make clean
	cd src/samtools && make clean
	cd src/bcftools && make clean
	cd src/bamStats && make clean
	mv src/homer/configureHomer.pl . && rm -rf src/homer/ && mkdir -p src/homer/ && mv configureHomer.pl src/homer/
	rm -rf $(TARGETS) $(TARGETS:=.o) src/FastQC src/cutadapt/ src/venv/ src/java/ src/weblogo/ src/gs/ src/blat/
