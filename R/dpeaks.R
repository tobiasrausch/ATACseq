library(DESeq2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(genefilter)

# Load count matrix and sample info
args=commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
colnames(x) = gsub(".final", "", colnames(x))
peaks = x[,1:4]
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

# Select ntop most variable peaks
ntop = 2500
rv = rowVars(assay(vsd))
selectIDX = order(rv, decreasing=T)[seq_len(min(ntop, length(rv)))]

# Sample heatmap
sampleDists = dist(t(assay(vsd)[selectIDX,]))
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(vsd$name)
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col=colors)

# PCA
pca = prcomp(t(assay(vsd)[selectIDX,]))
print(summary(pca))

# Check loadings (which peaks contribute most to each PC)
loadings = abs(pca$rotation)
contribution = as.data.frame(sweep(loadings, 2, colSums(loadings), "/"))
contribution = contribution[with(contribution, order(-PC1)),]
print(head(contribution))

# Plot PCA
pcaData = as.data.frame(pca$x)
pcaData$name=rownames(pcaData)
pcaData=merge(pcaData, s)
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))
p=ggplot(data=pcaData, aes(x = PC1, y = PC2, color=condition)) + geom_point(size=3)
p=p+xlab(paste0("PC1: ", percentVar[1], "% variance"))
p=p+ylab(paste0("PC2: ", percentVar[2], "% variance"))
p
q=ggplot(data=pcaData, aes(x = PC3, y = PC4, color=condition)) + geom_point(size=3)
q=q+xlab(paste0("PC3: ", percentVar[3], "% variance"))
q=q+ylab(paste0("PC4: ", percentVar[4], "% variance"))
q

# Differential expression
dds = DESeq(dds)

# log-foldchange > 1, adjusted p-value < 5%
#res = results(dds, lfcThreshold=0.25, alpha=0.05)
res = results(dds, lfcThreshold=0, alpha=0.05)
print(mcols(res, use.names=T))
print(summary(res))

# Histogram
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white", xlim=c(0,1), main="Histogram of p-values", xlab="p-value")

# MA-plot
plotMA(res, ylim = c(-5, 5))

# Significant results
resSig = subset(res, padj<0.05)
resSig = resSig[order(resSig$pvalue),]

# Heatmap
mat = assay(dds)[rownames(resSig),]
mat = mat - rowMeans(mat)
anno = as.data.frame(colData(dds)[, c("name", "condition")])
rownames(mat) = NULL
pheatmap(mat, annotation_col = anno, scale="row")      

# Write results
resSig = as.data.frame(resSig)
resSig$id = rownames(resSig)
resSig = merge(resSig, peaks)
write.table(resSig, file="res_sig.tsv", col.names=T, row.names=F, quote=F, sep="\t")
