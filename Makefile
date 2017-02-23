# External sources
PBASE=$(shell pwd)

# Targets
TARGETS = .fastqc .weblogo .blat .homer .picard 

all:   	$(TARGETS)

.fastqc:
	cd src && wget 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip' && unzip fastqc_v0.11.5.zip && chmod 755 FastQC/fastqc && rm fastqc_v0.11.5.zip && cd ../ && touch .fastqc

.weblogo:
	cd src/ && wget 'http://weblogo.berkeley.edu/release/weblogo.2.8.2.tar.gz' && tar -xzf weblogo.2.8.2.tar.gz && rm -rf weblogo.2.8.2.tar.gz && cd ../ && touch .weblogo

.blat: .blat
	mkdir src/blat/ && cd src/blat/ && wget 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat' && chmod 755 blat && cd ../../ && touch .blat

.homer: .weblogo .blat
	export PATH=${PBASE}/src/gs/bin:${PBASE}/src/weblogo:${PBASE}/src/blat:${PATH} && cd src/homer/ && perl configureHomer.pl -install homer && perl configureHomer.pl -install hg19 && cd ../../ && touch .homer

.picard: $(PICARDSOURCES)
	cd src/picard && ./gradlew shadowJar && cd ../../ && touch .picard

clean:
	cd src/picard && ./gradlew clean
	mv src/homer/configureHomer.pl . && rm -rf src/homer/ && mkdir -p src/homer/ && mv configureHomer.pl src/homer/
	rm -rf $(TARGETS) $(TARGETS:=.o) src/FastQC src/weblogo/ src/blat/
