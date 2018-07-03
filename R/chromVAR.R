library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BiocParallel)
library(DESeq2)


register(MulticoreParam(8))

#bamfiles = Sys.glob(c("./*U38*/*.bam"))
bamfiles = Sys.glob(c("PDXA5INI2AT100517/PDXA5INI2AT100517.final.bam", "PDXCoALL8INIAT190118/PDXCoALL8INIAT190118.final.bam", "PDXU38INIAT190118/PDXU38INIAT190118.final.bam", "PDXU58INIAT261017/PDXU58INIAT261017.final.bam"))

peaks = getPeaks("peaks.bed", sort_peaks=TRUE)

# Resize peaks
peaks = resize(peaks, width=500, fix="center")

# Fragment counting and peak filtering
fragment_counts = getCounts(bamfiles, peaks, paired=TRUE, by_rg = FALSE, format = "bam", colData = DataFrame(sample=c("A5", "CoALL8", "U38", "U58")))
counts_gc = addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg19)
counts_filt = filterSamples(counts_gc, min_depth=1500, min_in_peaks=0.05, shiny=FALSE)
counts_filtered = filterPeaks(counts_filt, non_overlapping=TRUE)
peakIdx = filterPeaks(counts_filt, non_overlapping=TRUE, ix_return=TRUE)
peaks = peaks[peakIdx]

# Motif search
motifs = getJasparMotifs()
motif_ix = matchMotifs(motifs, counts_filtered, genome=BSgenome.Hsapiens.UCSC.hg19, out="scores")

# Summary info
print(motif_ix)
print(dim(assays(motif_ix)$motifMatches))
print(dim(assays(motif_ix)$motifScores))
print(dim(assays(motif_ix)$motifCounts))

# Motif positions
motif_pos = matchMotifs(motifs, counts_filtered, genome=BSgenome.Hsapiens.UCSC.hg19, out="positions")

# Motif enrichment
bg = getBackgroundPeaks(object = counts_filtered)
dev = computeDeviations(object = counts_filtered, annotations = motif_ix, background_peaks = bg)

# Variability
variability = computeVariability(dev)
plotVariability(variability, n=50, use_plotly = FALSE)
write.table(variability, file="variability.tsv", col.names=T, row.names=F, quote=F, sep="\t")

# Inspect dev
tfs = unique(rowData(dev)$name)
print(head(tfs))

# t-SNE plots
tsne_results = deviationsTsne(dev, threshold = 1.5, perplexity = 0.9)
tsne_plots = plotDeviationsTsne(dev, tsne_results, annotation_name = "CTCF", sample_column = "sample", shiny=F)
tsne_plots[[1]]
tsne_plots[[2]]


# Plot TFs
tsne_plots = plotDeviationsTsne(dev, tsne_results, annotation_name = "GATA1::TAL1", sample_column = "sample", shiny=F)
tsne_plots[[2]]
tsne_plots = plotDeviationsTsne(dev, tsne_results, annotation_name = "NFIC::TLX1", sample_column = "sample", shiny=F)
tsne_plots[[2]]
tsne_plots = plotDeviationsTsne(dev, tsne_results, annotation_name = "HOXA5", sample_column = "sample", shiny=F)
tsne_plots[[2]]
tsne_plots = plotDeviationsTsne(dev, tsne_results, annotation_name = "NKX2-3", sample_column = "sample", shiny=F)
tsne_plots[[2]]
