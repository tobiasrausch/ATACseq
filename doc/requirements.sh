#!/bin/bash

grep -P "library\([^)]*\)" *.Rmd | sed 's/^.*library(\([^)]*\).*$/"\1"/'  | sort | uniq | tr '\n' ',' | sed 's/^/install.packages("BiocManager", repos="http:\/\/cran.us.r-project.org");BiocManager::install(c("bookdown",/' | sed 's/,$/), ask=F);/' > requirements.R
