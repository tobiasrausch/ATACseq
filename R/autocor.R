library(ggplot2)
library(scales)
library(reshape2)

args=commandArgs(trailingOnly=TRUE)
ac=read.table(args[1], skip=1, header=F)
colnames(ac)=c("dist", "Same", "Opposite")
ac=melt(ac, id.vars=c("dist"))

# Plot
pdf(paste0(args[1], ".pdf"))
p1=ggplot(data=ac, aes(x=dist, y=value))
p1=p1 + geom_line(aes(color=variable))
p1=p1 + xlab("Relative Distance between Reads (bp)") + ylab("Total Read Pairs") + ggtitle("Autocorrelation Analysis")
p1=p1 + labs(color="Strand")
p1
dev.off()
print(warnings())

