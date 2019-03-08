library(ggplot2)
library(scales)
library(corrplot)

args=commandArgs(trailingOnly=TRUE)
x=read.table(args[1], header=T)

colnames = c("BpCov1ToCov2Ratio","BpCov1ToCovNRatio","DuplicateFraction","ErrorRate","FRiP","FilteredPeaks","FractionChrM","FractionPeaksRetained","MappedSameChrFraction","PeakSaturation","SDCoverage","UnfilteredPeaks","UnmappedFraction","MappedReads", "TssEnrichment")

if (!("Study" %in% colnames(x))) { x$group = x$Sample; x$Study = "Samples"} else { x$group = x$Study; }

# Boxplots
for (col in colnames) {
 p = ggplot(data=x, aes_string(x="group", y=col)) + geom_boxplot(aes(group=group)) + scale_y_continuous(labels=comma) + coord_flip()
 ggsave(paste0(col, ".boxplot.png"))
}

# Correlation plot
for (s in x$Study) {
 df = x[x$Study==s,]
 df = df[, !(names(df) %in% c("Sample","Study", "group"))]

 # Cor-Plot 
 cr = cor(df)
 png(paste0(s, ".png"))
 corrplot(cr, type="upper", order="hclust")
 dev.off()
}
