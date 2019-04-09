library(karyoploteR)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=F)
x = x[,1:3]
colnames(x) = c("chrom", "start","end")
gr = toGRanges(x)

png("karyoplot.png")
kp = plotKaryotype(genome="hg19") #<- use the genome you need, if not human
kp = kpPlotDensity(kp, data=gr, col="blue")
dev.off()
