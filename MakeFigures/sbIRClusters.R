library(rtracklayer)
library(Rsamtools)
library(grid)
library(GenomicAlignments)
library(ggplot2)
library(GGally)
library(edgeR)
library(stringr)
library(EDASeq)
library(dplyr)
library(matrixStats)
library(gridExtra)
library(reshape2)
library(scales)
library(bigPint)

source("functions.R")

data("soybean_ir")
data <- soybean_ir

# Filter, normalize, and standardize the data so each gene has mean=0 and stdev=1
res <- filterStandardize(data)
# Fitered data standardized
filts <- res[["filts"]]
# Non-filtered data standardized
datas <- res[["datas"]]
# Hierarchical clustering object
hc <- res[["hc"]]
# Full data standardized
fulls <- rbind(datas, filts)

# Number of clusters
nC = 4
# Number of samples
nCol = 6
colList = scales::hue_pal()(nC+1)
colList <- colList[c(4, 3, 2, 5, 1)]
# Hierarchical clustering
k = cutree(hc, k=nC)

yMin = min(datas[,1:nCol])
yMax = max(datas[,1:nCol])

# Create background boxplot data
boxDat <- melt(fulls, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")

# Plot parallel coordinate lines over boxplots for non-filtered data
plot_clusters = lapply(1:nC, function(i){
  x = as.data.frame(datas[which(k==i),])
  nGenes = nrow(x)
  x$cluster = "color"
  x$cluster2 = factor(x$cluster)
  xNames = rownames(x)

  x$ID = xNames
  pcpDat <- melt(x[,c(1:(nCol+1))], id.vars="ID")
  colnames(pcpDat) <- c("ID", "Sample", "Count")
  boxDat$Sample <- as.character(boxDat$Sample)
  pcpDat$Sample <- as.character(pcpDat$Sample)
  
  p = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[i+1], alpha=0.05) + ylab("Standardized Count") + ggtitle(paste("Cluster ", i, " Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=18), axis.text=element_text(size=18), axis.title=element_text(size=18))
  p
})
# Create the filtered cluster plot
filts = as.data.frame(filts)
filts$cluster = "color"
filts$cluster2 = factor(filts$cluster)
nGenes = nrow(filts)
filts$ID = rownames(filts)
colnames(filts)[1:6] = colnames(datas)[1:6]
pcpDat <- melt(filts[,c(1:7)], id.vars="ID")
colnames(pcpDat) <- c("ID", "Sample", "Count")
pcpDat$Sample <- as.character(pcpDat$Sample)
plot_filtered = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[1], alpha=0.05) + ylab("Standardized Count") + ggtitle(paste("Filtered Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=18), axis.text=element_text(size=18), axis.title=element_text(size=18))

# Plot all clusters as a grid
do.call("grid.arrange", c(append(plot_clusters, list(plot_filtered)), ncol=ceiling(nC/2)))
