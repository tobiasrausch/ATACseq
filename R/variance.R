library(pvca)
library(DESeq2)
library(variancePartition)

x = read.table("sigcount.notss.counts.gz", header=T)
s = read.table("sample.info.tsv", header=T)
rownames(x) = x$id
rownames(s) = s$name

dds = DESeqDataSetFromMatrix(countData = x[,5:ncol(x)], colData = s, design = ~ 1)
dds = estimateSizeFactors(dds)
form = ~ (1|celltype) + (1|sex)
vsd = vst(dds)
varPart = fitExtractVarPartModel( assay(vsd), form, s)
vp = sortCols(varPart)
plotPercentBars(vp[1:10,])
plotVarPart(vp)


# pvca
pct_threshold = 0.5
batch.factors = c("celltype", "sex")
#batch.factors = c("patient", "sex")
es = ExpressionSet(assay(vsd), phenoData=AnnotatedDataFrame(s), featureData=AnnotatedDataFrame(x[,1:4]))
pvcaObj = pvcaBatchAssess(es, batch.factors, pct_threshold)

png("weightedVariance.png")
bp = barplot(pvcaObj$dat, xlab = "Effects", ylab = "Weighted average proportion variance", ylim = c(0, 1.1), col=c("blue"), las=2, main="PVCA estimation bar chart")
axis(1, at=bp, labels=pvcaObj$label, xlab="Effects", cex.axis=0.5, las=2)
values = pvcaObj$dat
new_values = round(values, 3)
text(bp, pvcaObj$dat, labels=new_values, pos=3, cex=0.8)
dev.off()
