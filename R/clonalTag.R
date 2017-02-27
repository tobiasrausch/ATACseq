library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
ct=read.table(args[1], skip=1, header=F)

# Plot
pdf(paste0(args[1], ".pdf"))
p1=ggplot(data=ct, aes(x=V1, y=V2))
p1=p1 + geom_bar(stat="identity", color="lightblue")
p1=p1 + xlab("Reads per position") + ylab("Fraction of total Reads") + ggtitle("Clonal Tag Distribution")
p1
dev.off()
print(warnings())

