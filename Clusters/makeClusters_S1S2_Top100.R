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

# Read in data
data(soybean_cn)
data <- soybean_cn
data(soybean_cn_metrics)
metrics <- soybean_cn_metrics

# Only focus on S1_S2
data <- data[,c(1:7)]
colnames(data) <- c("ID","S1.1","S1.2","S1.3","S2.1","S2.2","S2.3")
metrics <- metrics["S1_S2"]

boxDat <- data %>% gather(key, val, c(-ID))
colnames(boxDat) <- c("ID", "Sample", "Count")

baseOutDir = "Clustering_S1_S2_Top100"

id100 = metrics[["S1_S2"]]$ID[1:100]
id100df = which(data$ID %in% id100)
data = data[id100df,]
#data = data[,-1]

plotName = "S1_S2"
outDir = "Clustering_S1_S2"
outDir= paste0(getwd(),"/",baseOutDir)

dendo = data[,-1]
rownames(data[,-1]) = NULL
# Euclidean distance between rows of matrix
d = dist(as.matrix(dendo))
# Hierarchical clustering using ward.D linkage
hc = hclust(d, method="ward.D")
plotName = "S1_S2_Dendogram"
jpeg(file = paste(outDir, "/", plotName, ".jpg", sep=""))
plot(hc,main=plotName, xlab=NA, sub=NA)
invisible(dev.off())

getPCP <- function(nC){
  
  set.seed(1)
  colList = c("red","magenta", "darkorange","darkgreen","blue","purple")
  k = cutree(hc, k=nC)
  ###########################
  plot_clusters = lapply(1:nC, function(i){
    plotName = "S1_S2"
    x = as.data.frame(data[which(k==i),])
    nGenes = nrow(x)
    x$cluster = "color"
    x$cluster2 = factor(x$cluster)
    xNames = rownames(x)
    write.table(xNames, file = paste(outDir, "/", plotName, "_", nC, "_", i, ".txt", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)  
    
    xPCP <- x[,c(1:7)]
    pcpDat2 <- xPCP %>% gather(key, val,c(-ID))
    colnames(pcpDat2) <- c("ID", "Sample", "Count")
    p <- ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat2, aes_string(x = 'Sample', y = 'Count', group = 'ID'), color = colList[i]) + xlab(paste("Cluster ", i, " (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(legend.position = "none", axis.title=element_text(size=12), axis.text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=15), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))
    fileName = paste(outDir, "/", plotName, "_", nC, "_", i, ".jpg", sep="")
    jpeg(fileName)
    plot(p)
    invisible(dev.off())
    p    
  })
  ###########################
  jpeg(file = paste(outDir, "/S1_S2_", nC, ".jpg", sep=""), height = 200 * ceiling((nC+1)/3), width = min(200 * (nC+1), 800))
  # We allow up to 4 (now 3) plots in each column
  p = do.call("grid.arrange", c(plot_clusters, ncol=3)) #change from 4 to 3 ceiling(nC/3)
  invisible(dev.off())
}

for (i in c(2,3,4,5,6)){
  getPCP(i)
}
