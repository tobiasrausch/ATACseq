ATAC = atac.html atac.R
TARGETS = $(ATAC)

PBASE = $(shell pwd)

all: $(TARGETS)

atac: $(ATAC)

%.html: %.Rmd
	export PATH=${PBASE}/R/bin:${PATH} && Rscript -e 'rmarkdown::render("$<")'

%.R: %.Rmd
	export PATH=${PBASE}/R/bin:${PATH} && Rscript -e 'knitr::purl("$<")'

book: *.Rmd
	export PATH=${PBASE}/R/bin:${PATH} && Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::gitbook")'

.Rinstall:
	wget 'http://cran.rstudio.com/src/base/R-latest.tar.gz' && tar -xzf R-latest.tar.gz && cd R-* && ./configure --prefix=$(PBASE)/R && make && make install && cd ../ && rm -rf R-latest.tar.gz R-*/ && touch .Rinstall

requirements.R: .Rinstall
	export PATH=${PBASE}/R/bin:${PATH} && ./requirements.sh

install: .Rinstall requirements.R
	export PATH=${PBASE}/R/bin:${PATH} && Rscript requirements.R

clean:
	rm -rf $(TARGETS) Rplots.pdf atac_files/ _book/ _bookdown_files/ requirements.R _main.rds

distclean:
	rm -rf .Rinstall R/

.PHONY: distclean clean all
