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

baseOutDir = "Clustering_S1_S2"

rowNames <- data[,1]
data <- data[,-1]
data <- DGEList(counts=data)

# Filter
data <- data[rowSums(data$counts>1)>=ncol(data)/2,]
cpm.data.new <- cpm(data, TRUE, TRUE)
cpm.data.norm <- betweenLaneNormalization(cpm.data.new, which="full", round=FALSE)
data = cpm.data.norm

RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

datat = data
data = as.data.frame(datat)
data = mutate(data, mean = (S1.1+S1.2+S1.3+S2.1+S2.2+S2.3)/ncol(data), stdev = RowSD(cbind(S1.1,S1.2,S1.3,S2.1,S2.2,S2.3)))
rownames(data)=rownames(datat)

q1T = as.numeric(summary(data$mean)["1st Qu."])
dataq1 = subset(data,mean>q1T)
q1Ts = as.numeric(summary(dataq1$stdev)["1st Qu."])
dataq1 = subset(dataq1,stdev>q1Ts)
filt = subset(data,mean<=q1T|stdev<=q1Ts)
model = loess(mean ~ stdev, data=dataq1)
dataq1 = dataq1[which(sign(model$residuals) == 1),]
dataq1 = dataq1[,1:(ncol(dataq1)-2)]
dataq1s = t(apply(as.matrix(dataq1), 1, scale))
colnames(dataq1s)=colnames(dataq1)
colnames(dataq1)=colnames(dataq1) 
filt = filt[,1:(ncol(filt)-2)]
colnames(filt)=colnames(dataq1)
filt = rbind(filt,dataq1[which(sign(model$residuals) == -1),])
filts = t(apply(as.matrix(filt), 1, scale))
colnames(filts)=colnames(dataq1)
colnames(filt)=colnames(dataq1)

plotName = "S1_S2"
outDir = "Clustering_S1_S2"
outDir= paste0(getwd(),"/",outDir)

dendo = dataq1s
rownames(dataq1s) = NULL
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
  colList = scales::hue_pal()(nC+1)
  k = cutree(hc, k=nC)
  ###########################
  plot_clusters = lapply(1:nC, function(i){
    plotName = "S1_S2"
    x = as.data.frame(dataq1s[which(k==i),])
    nGenes = nrow(x)
    x$cluster = "color"
    x$cluster2 = factor(x$cluster)
    xNames = rownames(x)
    write.table(xNames, file = paste(outDir, "/", plotName, "_", nC, "_", i, ".txt", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)  
    
    p = ggparcoord(x, columns=1:6, groupColumn=8, scale="globalminmax", alphaLines = 0.2) + xlab(paste("Cluster ", i, " (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count") + theme(legend.position = "none", axis.title=element_text(size=12), axis.text=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = c("color" = colList[i+1]))
    fileName = paste(outDir, "/", plotName, "_", nC, "_", i, ".jpg", sep="")
    jpeg(fileName)
    plot(p)
    invisible(dev.off())
    p
  })
  ###########################
  filts = as.data.frame(filts)
  filts$cluster = "color"
  filts$cluster2 = factor(filts$cluster)
  nGenes = nrow(filts)
  plot_filtered = ggparcoord(filts, columns=1:6, groupColumn=8, scale="globalminmax", alphaLines = 0.01) + xlab(paste("Filtered (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count") + theme(legend.position = "none", axis.title=element_text(size=12), axis.text=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = c("color" = colList[1]))
  ###########################
  jpeg(file = paste(outDir, "/S1_S2_", nC, ".jpg", sep=""), height = 200 * ceiling((nC+1)/3), width = min(200 * (nC+1), 800))
  # We allow up to 4 (now 3) plots in each column
  p = do.call("grid.arrange", c(append(plot_clusters, list(plot_filtered)), ncol=3)) #change from 4 to 3 ceiling(nC/3)
  invisible(dev.off())
}

for (i in c(2,3,11)){
  getPCP(i)
}
