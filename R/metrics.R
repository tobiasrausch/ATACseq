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

# Correlation plot
for (s in x$Study) {
 df = x[x$Study==s,]
 df = df[, !(names(df) %in% c("Sample","Study"))]

 # Cor-Plot 
 cr = cor(df)
 png(paste0(s, ".png"))
 corrplot(cr, type="upper", order="hclust")
 dev.off()
}
quit()

# Correlation between any 2 QC metrics
for (s in x$Study) {
 df = x[x$Study==s,]
 for(c1 in colnames) {
   for(c2 in colnames) {
     if (c1 != c2) {
      dfsub = df[,names(df) %in% c(c1, c2)]
      png(paste0(s, ".", c1, ".", c2, ".png"))
      chart.Correlation(dfsub)
      dev.off()
    }
   }
 }
}
    
