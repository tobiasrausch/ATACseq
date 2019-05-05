#!/bin/bash

wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74912/suppl/GSE74912_ATACseq_All_Counts.txt.gz'

# Extract blood samples
cat sample.mapping | grep "5483$" | grep -v "CLP" > select.samples
cat select.samples | sed 's/Donor5483-//' | sed 's/-/_/' | sed 's/-[^\t]*//' > blood.samples

# Extract respective columns
COLS=`zcat GSE74912_ATACseq_All_Counts.txt.gz | head -n 1 | tr '\t' '\n' | grep -n -w -Ff <(cut -f 1 select.samples) | sed 's/:.*//' | tr '\n' ',' | sed 's/,$//'`
zcat GSE74912_ATACseq_All_Counts.txt.gz | cut -f 1-3,${COLS} | sed 's/Donor5483-//g' | sed 's/-\([^-]*\)-[^\t]*/_\1/g' | gzip -c > atac.data.gz
rm select.samples
