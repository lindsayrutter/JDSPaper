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

data("soybean_ir")
data("soybean_ir_metrics")
data <- soybean_ir
metrics <- soybean_ir_metrics[["N_P"]]

logSoy = soybean_ir
logSoy[,-1] <- log(soybean_ir[,-1]+1)

# Make sure each gene has at least one count in at least half of the six samples
filterLow = which(rowSums(data[,-1])<=ncol(data[,-1])/2)
filt1 <- data[filterLow,]
data <- data[-filterLow,]

RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

rownames_filt1 <- filt1$ID
filt1 <- filt1[,-1]
filt1 = mutate(filt1, mean = (N.1+N.2+N.3+P.1+P.2+P.3)/6, stdev = RowSD(cbind(N.1,N.2,N.3,P.1,P.2,P.3)))
rownames(filt1) <- rownames_filt1

data_Rownames <- data$ID
data = data[,-1]
rownames(data) <- data_Rownames

cpm.data.new <- cpm(data, TRUE, TRUE)
# Between-lane normalization for sequencing depth and possibly other distributional differences between lanes.
cpm.data.norm <- betweenLaneNormalization(cpm.data.new, which="full", round=FALSE)
data = cpm.data.norm

data = as.data.frame(data)
# Mutate is to add columns to a data frame
data = mutate(data, mean = (N.1+N.2+N.3+P.1+P.2+P.3)/6, stdev = RowSD(cbind(N.1,N.2,N.3,P.1,P.2,P.3)))
rownames(data)=data_Rownames

# The first quartile threshold of mean counts across the 5 samples
q1T = as.numeric(summary(data$mean)["1st Qu."])
dataq1 = subset(data,mean>q1T)
q1Ts = as.numeric(summary(dataq1$stdev)["1st Qu."])
dataq1 = subset(dataq1,stdev>q1Ts)
filt = subset(data,mean<=q1T|stdev<=q1Ts)
filt = rbind(filt, filt1)

model = loess(mean ~ stdev, data=dataq1)
dataq1p = dataq1[which(sign(model$residuals) == 1),]
dataq1n = dataq1[which(sign(model$residuals) == -1),]
dataq1p = dataq1p[,1:6]
# Scale each row, so now negative and positive values
dataq1s = t(apply(as.matrix(dataq1p), 1, scale))
colnames(dataq1s)=colnames(dataq1p)
filt = filt[,1:6]
dataq1n = dataq1n[,1:6]
# filt is 31394
filt = rbind(filt,dataq1n)
filts = t(apply(as.matrix(filt), 1, scale))

dendo = dataq1s
rownames(dendo) = NULL
d = dist(as.matrix(dendo))
hc = hclust(d, method="ward.D")

plotName = "N_P"
outDir = "Clustering_data"

fileName = paste(getwd(), "/", outDir, "/dendodgram.jpg", sep="")
jpeg(fileName)
plot(hc,main="data Dendogram", xlab=NA, sub=NA)
invisible(dev.off())

#####################################################

getPCP <- function(nC){

set.seed(1)
colList = scales::hue_pal()(nC+1)
k = cutree(hc, k=nC)

yMin = min(dataq1s[,1:6])
yMax = max(dataq1s[,1:6])

###########################

sbsDF <- data.frame()
for (i in 1:nC){
  x = as.data.frame(dataq1s[which(k==i),])
  xNames = rownames(x)
  xPValues = metrics[which(metrics$ID %in% xNames),]$PValue
  sbsDF = rbind(sbsDF, data.frame(Cluster = paste("Cluster", i), PValue = xPValues))
}

plot_clusters = lapply(1:nC, function(i){
  x = as.data.frame(dataq1s[which(k==i),])
  nGenes = nrow(x)
  x$cluster = "color"
  x$cluster2 = factor(x$cluster)
  xNames = rownames(x)
  xPValues = metrics[which(metrics$ID %in% xNames),]$PValue
  scatMatMetrics = list()
  scatMatMetrics[["N_P"]] = metrics[which(metrics$ID %in% xNames),]
  scatMatMetrics[["N_P"]]$PValue = 10e-10
  scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
  
  fileName = paste(getwd(), "/", outDir, "/", "SM_", nC, "_", i, ".jpg", sep="")
  plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "PValue", threshVal = 0.05/nrow(logSoy), degPointColor = colList[i+1], fileName=fileName)

  p = ggparcoord(x, columns=1:6, groupColumn=8, scale="globalminmax", alphaLines = 0.1) + xlab(paste("Cluster ", i, " (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count") + theme(legend.position = "none", axis.title=element_text(size=11), axis.text=element_text(size=11), axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = c("color" = colList[i+1])) + ylim(yMin, yMax)
  fileName = paste(getwd(), "/", outDir, "/", plotName, "_", nC, "_", i, ".jpg", sep="")
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
  
  xNames = rownames(filts)
  xPValues = metrics[which(metrics$ID %in% xNames),]$PValue
  sbsDF = rbind(sbsDF, data.frame(Cluster = paste("Filtered"), PValue = xPValues))
  
  ggBP = ggplot(sbsDF, aes(x=Cluster, y=PValue)) +
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(outlier.shape=NA, aes(fill=Cluster), alpha = 0.3) +
  geom_point(aes(fill=Cluster), shape=21, position=position_jitter(width=0.3), alpha=0.1) +
  scale_fill_manual(values=colList[c(2:length(colList), 1)])
  jpeg(file = paste(getwd(), "/", outDir, "/boxplot_", nC, ".jpg", sep=""), width=1000, height=700)
  ggBP
  invisible(dev.off())
  
  scatMatMetrics = list()
  scatMatMetrics[["N_P"]] = metrics[which(metrics$ID %in% xNames),]
  scatMatMetrics[["N_P"]]$PValue = 10e-10
  scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
  fileName = paste(getwd(), "/", outDir, "/", "SM_", nC, "_filtered.jpg", sep="")
  
  plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "PValue", threshVal = 0.05/nrow(logSoy), degPointColor = colList[1], fileName=fileName)
  
  plot_filtered = ggparcoord(filts, columns=1:6, groupColumn=8, scale="globalminmax", alphaLines = 0.01) + xlab(paste("Filtered (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count") + theme(legend.position = "none", axis.title=element_text(size=11), axis.text=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = c("color" = colList[1])) + ylim(yMin, yMax)

  jpeg(file = paste(getwd(), "/", outDir, "/", plotName, "_", nC, ".jpg", sep=""), width=1000, height=700)
  # We allow up to 4 plots in each column
  p = do.call("grid.arrange", c(append(plot_clusters, list(plot_filtered)), ncol=ceiling(nC/2)))
  invisible(dev.off())
  
plot_clustersSig = lapply(1:nC, function(i){ 
  x = as.data.frame(dataq1s[which(k==i),])
  x$cluster = "color"
  x$cluster2 = factor(x$cluster)
  xNames = rownames(x)
  metricPValue = metrics[which(as.character(metrics$ID) %in% xNames),]
  sigID = metricPValue[metricPValue$PValue<0.05/nrow(soybean_ir),]$ID
  xSig = x[which(rownames(x) %in% sigID),]
  xSigNames = rownames(xSig)
  saveRDS(xSigNames, file=paste0(getwd(), "/", outDir, "/Sig_", nC, "_", i, ".Rds"))

  if (nrow(xSig)>0){
    scatMatMetrics = list()
    scatMatMetrics[["N_P"]] = metrics[which(metrics$ID %in% xSigNames),]
    scatMatMetrics[["N_P"]]$PValue = 10e-10
    scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
    fileName = paste(getwd(), "/", outDir, "/", "SM_Sig_", nC, "_", i, ".jpg", sep="")
    plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "PValue", threshVal = 0.05/nrow(logSoy), degPointColor = colList[i+1], fileName=fileName)
    pSig = ggparcoord(xSig, columns=1:6, groupColumn=8, scale="globalminmax", alphaLines = 1) + xlab(paste("Cluster ", i, " (n=", format(nrow(xSig), big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count") + theme(legend.position = "none", axis.title=element_text(size=11), axis.text=element_text(size=11), axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = c("color" = colList[i+1])) + ylim(yMin, yMax)
  }else{
    scatMatMetrics = list()
    scatMatMetrics[["N_P"]] = metrics[1,]
    scatMatMetrics[["N_P"]]$PValue = 1
    scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
    fileName = paste(getwd(), "/", outDir, "/", "SM_Sig_", nC, "_", i, ".jpg", sep="")
    plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "PValue", threshVal = 0.05/nrow(logSoy), degPointColor = colList[i+1], fileName=fileName)
    xSig = x[1,]
    pSig = ggparcoord(xSig, columns=1:6, groupColumn=8, scale="globalminmax", alphaLines = 0) + xlab(paste("Cluster ", i, " (n=", format(nrow(xSig)-1, big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count") + theme(legend.position = "none", axis.title=element_text(size=11), axis.text=element_text(size=11), axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = c("color" = colList[i+1])) + ylim(yMin, yMax)
  }
  fileName = paste(getwd(), "/", outDir, "/", plotName, "_Sig_", nC, "_", i, ".jpg", sep="")
  jpeg(fileName)
  plot(pSig)
  invisible(dev.off())
  pSig
})
  
  xNames = rownames(filts)
  metricPValue = metrics[which(as.character(metrics$ID) %in% xNames),]
  sigID = metricPValue[metricPValue$PValue<0.05/nrow(soybean_ir),]$ID
  filtsSig = filts[which(rownames(filts) %in% sigID),]
  filtsSigNames = rownames(filtsSig)
  saveRDS(filtsSigNames, file=paste0(getwd(), "/", outDir, "/Sig_", nC, "_Filtered.Rds"))
  
  scatMatMetrics = list()
  scatMatMetrics[["N_P"]] = metrics[which(metrics$ID %in% filtsSigNames),]
  scatMatMetrics[["N_P"]]$PValue = 10e-10
  scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
  fileName = paste(getwd(), "/", outDir, "/", "SM_Sig_", nC, "_filtered.jpg", sep="")
  plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "PValue", threshVal = 0.05/nrow(logSoy), degPointColor = colList[1], fileName=fileName)
  
  if (nrow(filtsSig)>0){
    plot_filteredSig = ggparcoord(filtsSig, columns=1:6, groupColumn=8, scale="globalminmax", alphaLines = 1) + xlab(paste("Filtered(n=", format(nrow(filtsSig), big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count") + theme(legend.position = "none", axis.title=element_text(size=11), axis.text=element_text(size=11), axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = c("color" = colList[1])) + ylim(yMin, yMax)
  } else{
    pSig = filtsSig[1,]
    plot_filteredSig = ggparcoord(filtsSig, columns=1:6, groupColumn=8, scale="globalminmax", alphaLines = 0) + xlab(paste("Filtered (n=", format(nrow(filtsSig), big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count") + theme(legend.position = "none", axis.title=element_text(size=11), axis.text=element_text(size=11), axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_manual(values = c("color" = colList[1])) + ylim(yMin, yMax)
  }

  jpeg(file = paste(getwd(), "/", outDir, "/", plotName, "_Sig_", nC, ".jpg", sep=""), width=1000, height=700)
  # We allow up to 4 plots in each column
  p = do.call("grid.arrange", c(append(plot_clustersSig, list(plot_filteredSig)), ncol=ceiling(nC/2)))
  invisible(dev.off())
}

for (i in 2:5){
  getPCP(i)
}
