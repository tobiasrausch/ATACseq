library(ggplot2)
library(scales)
library(reshape2)

args=commandArgs(trailingOnly=TRUE)
nf=read.table(args[1], header=T)
nf=nf[,1:5]
nf=melt(nf, id.vars=c("Offset"))

# Plot
pdf(paste0(args[1], ".pdf"))
p1=ggplot(data=nf, aes(x=Offset, y=value))
p1=p1 + geom_line(aes(color=variable))
p1=p1 + xlab("Distance from 5' end of Reads") + ylab("Nucleotide Frequency") + ggtitle("Genomic Nucleotide Frequency Relative to Read Position")
p1=p1 + labs(color="Nucleotide")
p1
dev.off()
print(warnings())

