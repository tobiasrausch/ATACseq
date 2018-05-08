library(ggplot2)
library(scales)
library(corrplot)
library(PerformanceAnalytics)
library(psych)

args=commandArgs(trailingOnly=TRUE)
x=read.table(args[1], header=T)

colnames = c("BpCov1ToCov2Ratio","BpCov1ToCovNRatio","DuplicateFraction","ErrorRate","FRiP","FilteredPeaks","FractionChrM","FractionPeaksRetained","MappedSameChrFraction","PeakSaturation","SDCoverage","UnfilteredPeaks","UnmappedFraction","MappedReads", "TssEnrichment")

if (!("Study" %in% colnames(x))) { x$Study = "Samples"; }

# Boxplots
for (col in colnames) {
 p = ggplot(data=x, aes_string(x="Study", y=col)) + geom_boxplot(aes(group=Study)) + scale_y_continuous(labels=comma)
 ggsave(paste0(col, ".boxplot.png"))
}

# Density curves
for (col in colnames) {
 p = ggplot(data=x, aes_string(x=col)) + geom_freqpoly(aes(group=Study, color=Study), bins=20) + scale_y_continuous(labels=comma) + scale_x_continuous(labels=comma)
 ggsave(paste0(col, ".density.png"))
}

quit()



df = x[x$Study=="Public",]
df = df[, !(names(df) %in% c("Sample","Study"))]

# Cor-Plot 
cr = cor(df)
corrplot(cr, type="upper", order="hclust")

# Correlation Chart
df = df[, !(names(df) %in% c("BpCov1ToCov2Ratio","BpCov1ToCovNRatio","DuplicateFraction","ErrorRate","FRiP","FilteredPeaks","FractionChrM","FractionPeaksRetained","MappedSameChrFraction","PeakSaturation","SDCoverage","Sample","UnfilteredPeaks","UnmappedFraction","Study"))]
chart.Correlation(df)
