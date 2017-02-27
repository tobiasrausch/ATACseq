library(ggplot2)
library(scales)
library(reshape2)

args=commandArgs(trailingOnly=TRUE)
genome=read.table(args[1], skip=1, header=F)
genome=genome[,c(1,3)]
colnames(genome)=c("gc","Genome")
gc=read.table(args[2], skip=1, header=F)
gc=gc[,c(1,3)]
colnames(gc)=c("gc", "ATACseq")
gc=merge(genome, gc)
gc=melt(gc, id.vars=c("gc"))

# Plot
pdf(paste0(args[2], ".pdf"))
p1=ggplot(data=gc, aes(x=gc, y=value))
p1=p1 + geom_line(aes(color=variable))
p1=p1 + xlab("GC-content of Fragments") + ylab("Normalized Fraction") + ggtitle("Fragment GC% Distribution")
p1=p1 + labs(color="Type")
p1
dev.off()
print(warnings())

