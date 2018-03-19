library(DESeq2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
		
# Load count matrix and sample info
args=commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
colnames(x) = gsub(".final", "", colnames(x))
s = read.table(args[2], header=T)
rownames(x) = x$id
rownames(s) = s$name
s$condition = factor(s$condition)

# Build DESeq object
dds = DESeqDataSetFromMatrix(countData = x[,5:ncol(x)], colData = s, design = ~ condition)
print(dds)

# Variance stabilizing transformation
vsd = vst(dds, blind=F)
meanSdPlot(assay(vsd), ranks=F)

# Sample heatmap
sampleDists = dist(t(assay(vsd)))
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(vsd$name)
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col=colors)

# PCA
pcaData = plotPCA(vsd, intgroup = c("condition"), returnData=T)
percentVar = round(100 * attr(pcaData, "percentVar"))
p=ggplot(data=pcaData, aes(x = PC1, y = PC2, color=condition)) + geom_point(size=3)
p=p+xlab(paste0("PC1: ", percentVar[1], "% variance"))
p=p+ylab(paste0("PC2: ", percentVar[2], "% variance"))
p

# Differential expression
dds = DESeq(dds)

# log-foldchange > 1, adjusted p-value < 5%
res = results(dds, lfcThreshold=1, alpha=0.05)
print(mcols(res, use.names=T))
print(summary(res))
