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

RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

# Make sure each gene has at least one count in at least half of the six samples
filterLow = which(rowSums(data[,-1])<=ncol(data[,-1])/2)
filt1 <- data[filterLow,]
rownames_filt1 <- filt1$ID
filt1 <- filt1[,-1]
filt1 = mutate(filt1, mean = (N.1+N.2+N.3+P.1+P.2+P.3)/6, stdev = RowSD(cbind(N.1,N.2,N.3,P.1,P.2,P.3)))
rownames(filt1) <- rownames_filt1

data <- data[-filterLow,]
data_Rownames <- data$ID
data = data[,-1]
rownames(data) <- data_Rownames
#Normalize and log
cpm.data.new <- cpm(data, TRUE, TRUE)
# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(cpm.data.new, which="full", round=FALSE)
data = as.data.frame(data)
# Add mean and standard deviation for each row/gene
data = mutate(data, mean = (N.1+N.2+N.3+P.1+P.2+P.3)/6, stdev = RowSD(cbind(N.1,N.2,N.3,P.1,P.2,P.3)))
rownames(data)=data_Rownames
data$ID <- data_Rownames
# Remove the genes with lowest quartile of mean and standard deviation
qT = as.numeric(summary(data$mean)["1st Qu."])
dataq = subset(data,mean>qT)
qTs = as.numeric(summary(dataq$stdev)["1st Qu."])
dataq = subset(dataq,stdev>qTs)
filt = subset(data,mean<=qT|stdev<=qTs)
filt <- rbind(filt[,-9], filt1)
filt$ID <- rownames(filt)

# Apply Loess model and further filter low gene counts
model = loess(mean ~ stdev, data=dataq)
dataqp = dataq[which(sign(model$residuals) == 1),]
dataqn = dataq[which(sign(model$residuals) == -1),]
dataqp = dataqp[,1:6]

#Scale filter data
filt = filt[,1:6]
filt = rbind(filt,dataqn[,1:6])

dataqps <- t(apply(as.matrix(dataqp[,1:6]), 1, scale))
filts <- t(apply(as.matrix(filt[,1:6]), 1, scale))
dataqps <- as.data.frame(dataqps)
colnames(dataqps) <- colnames(dataqp[,1:6])
dataqps$ID <- rownames(dataqps)
filts <- as.data.frame(filts)
colnames(filts) <- colnames(filt[,1:6])
filts$ID <- rownames(filts)
# Indices of the 9760 NAN rows. They had stdev=0 in the filt data
nID <- which(is.nan(filts$N.1))
# Set these filtered values that have all same values for samples to 0
filts[nID,1:6] <- 0

# Comine the filtered and remaining data
fulls <- rbind(dataqps, filts)
boxDat <- melt(fulls, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")
#ggplot(boxDat, aes(x=Sample, y=Count)) + geom_boxplot()

dendo = dataqps # or dataqps? (If do fulls, then have NAs introduced by conversion)
rownames(dendo) = NULL
d = dist(as.matrix(dendo))
hc = hclust(d, method="ward.D")

plotName = "N_P"
outDir = "Clustering_data_box3"

fileName = paste(getwd(), "/", outDir, "/dendodgram.jpg", sep="")
jpeg(fileName)
plot(hc, main="data Dendogram", xlab=NA, sub=NA)
invisible(dev.off())

logSoy = soybean_ir
logSoy[,-1] <- log(soybean_ir[,-1]+1)

#####################################################

getPCP <- function(nC){
  
  set.seed(1)
  colList = scales::hue_pal()(nC+1)
  k = cutree(hc, k=nC)
  
  yMin = min(dataqps[,1:6])
  yMax = max(dataqps[,1:6])
  
  ###########################
  
  sbsDF <- data.frame()
  for (i in 1:nC){
    x = as.data.frame(dataqps[which(k==i),])
    xNames = rownames(x)
    xFDR = metrics[which(metrics$ID %in% xNames),]$FDR
    sbsDF = rbind(sbsDF, data.frame(Cluster = paste("Cluster", i), FDR = xFDR))
  }
  
  plot_clusters = lapply(1:nC, function(i){
    x = as.data.frame(dataqps[which(k==i),])
    nGenes = nrow(x)
    x$cluster = "color"
    x$cluster2 = factor(x$cluster)
    xNames = rownames(x)
    xFDR = metrics[which(metrics$ID %in% xNames),]$FDR
    scatMatMetrics = list()
    scatMatMetrics[["N_P"]] = metrics[which(metrics$ID %in% xNames),]
    scatMatMetrics[["N_P"]]$FDR = 10e-10
    scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
    
    fileName = paste(getwd(), "/", outDir, "/", "SM_", nC, "_", i, ".jpg", sep="")
    plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "FDR", threshVal = 0.05, degPointColor = colList[i+1], fileName=fileName)
    
    x$ID = xNames
    
    pcpDat <- melt(x[,c(1:7)], id.vars="ID")
    colnames(pcpDat) <- c("ID", "Sample", "Count")
    boxDat$Sample <- as.character(boxDat$Sample)
    pcpDat$Sample <- as.character(pcpDat$Sample)

    p = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[i+1], alpha=0.2) + xlab(paste("Cluster ", i, " (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count")
    
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
  xFDR = metrics[which(metrics$ID %in% xNames),]$FDR
  sbsDF = rbind(sbsDF, data.frame(Cluster = paste("Filtered"), FDR = xFDR))
  
  ggBP = ggplot(sbsDF, aes(x=Cluster, y=FDR)) +
    stat_boxplot(geom ='errorbar') + 
    geom_boxplot(outlier.shape=NA, aes(fill=Cluster), alpha = 0.3) +
    geom_point(aes(fill=Cluster), shape=21, position=position_jitter(width=0.3), alpha=0.1) +
    scale_fill_manual(values=colList[c(2:length(colList), 1)])
  jpeg(file = paste(getwd(), "/", outDir, "/boxplot_", nC, ".jpg", sep=""), width=1000, height=700)
  ggBP
  invisible(dev.off())
  
  scatMatMetrics = list()
  scatMatMetrics[["N_P"]] = metrics[which(metrics$ID %in% xNames),]
  scatMatMetrics[["N_P"]]$FDR = 10e-10
  scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
  fileName = paste(getwd(), "/", outDir, "/", "SM_", nC, "_filtered.jpg", sep="")
  
  plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "FDR", threshVal = 0.05, degPointColor = colList[1], fileName=fileName)
  
  filts$ID = xNames
  colnames(filts)[1:6] = colnames(dataqps)[1:6]
  
  pcpDat <- melt(filts[,c(1:7)], id.vars="ID")
  colnames(pcpDat) <- c("ID", "Sample", "Count")
  pcpDat$Sample <- as.character(pcpDat$Sample)
  
  plot_filtered = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[1], alpha=0.2) + xlab(paste("Cluster ", i, " (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count")
  
  jpeg(file = paste(getwd(), "/", outDir, "/", plotName, "_", nC, ".jpg", sep=""), width=1000, height=700)
  # We allow up to 4 plots in each column
  p = do.call("grid.arrange", c(append(plot_clusters, list(plot_filtered)), ncol=ceiling(nC/2)))
  invisible(dev.off())
  
  plot_clustersSig = lapply(1:nC, function(i){ 
    x = as.data.frame(dataqps[which(k==i),])
    x$cluster = "color"
    x$cluster2 = factor(x$cluster)
    xNames = rownames(x)
    metricFDR = metrics[which(as.character(metrics$ID) %in% xNames),]
    sigID = metricFDR[metricFDR$FDR<0.05,]$ID
    xSig = x[which(rownames(x) %in% sigID),]
    xSigNames = rownames(xSig)
    nGenes = nrow(xSig)
    saveRDS(xSigNames, file=paste0(getwd(), "/", outDir, "/Sig_", nC, "_", i, ".Rds"))
    
    if (nrow(xSig)>0){
      scatMatMetrics = list()
      scatMatMetrics[["N_P"]] = metrics[which(metrics$ID %in% xSigNames),]
      scatMatMetrics[["N_P"]]$FDR = 10e-10
      scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
      fileName = paste(getwd(), "/", outDir, "/", "SM_Sig_", nC, "_", i, ".jpg", sep="")
      plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "FDR", threshVal = 0.05, degPointColor = colList[i+1], fileName=fileName)
      
      xSig$ID = xSigNames
      pcpDat <- melt(xSig[,c(1:7)], id.vars="ID")
      colnames(pcpDat) <- c("ID", "Sample", "Count")
      pcpDat$Sample <- as.character(pcpDat$Sample)
      
      pSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[i+1], alpha=0.2) + xlab(paste("Cluster ", i, " (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count")

    }else{
      scatMatMetrics = list()
      scatMatMetrics[["N_P"]] = metrics[1,]
      scatMatMetrics[["N_P"]]$FDR = 1
      scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
      fileName = paste(getwd(), "/", outDir, "/", "SM_Sig_", nC, "_", i, ".jpg", sep="")
      plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "FDR", threshVal = 0.05, degPointColor = colList[i+1], fileName=fileName)
      xSig = x[1,]

      xSig$ID = c("PlaceholderID")
      pcpDat <- melt(xSig[,c(1:7)], id.vars="ID")
      colnames(pcpDat) <- c("ID", "Sample", "Count")
      pcpDat$Sample <- as.character(pcpDat$Sample)
      
      pSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[i+1], alpha=0) + xlab(paste("Cluster ", i, " (n=0", ")",sep="")) + ylab("Count")
    }
    fileName = paste(getwd(), "/", outDir, "/", plotName, "_Sig_", nC, "_", i, ".jpg", sep="")
    jpeg(fileName)
    plot(pSig)
    invisible(dev.off())
    pSig
  })
  
  xNames = rownames(filts)
  metricFDR = metrics[which(as.character(metrics$ID) %in% xNames),]
  sigID = metricFDR[metricFDR$FDR<0.05,]$ID
  filtsSig = filts[which(rownames(filts) %in% sigID),]
  nGenes = nrow(filtsSig)
  filtsSigNames = rownames(filtsSig)
  saveRDS(filtsSigNames, file=paste0(getwd(), "/", outDir, "/Sig_", nC, "_Filtered.Rds"))
  
  scatMatMetrics = list()
  scatMatMetrics[["N_P"]] = metrics[which(metrics$ID %in% filtsSigNames),]
  scatMatMetrics[["N_P"]]$FDR = 10e-10
  scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
  fileName = paste(getwd(), "/", outDir, "/", "SM_Sig_", nC, "_filtered.jpg", sep="")
  plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "FDR", threshVal = 0.05, degPointColor = colList[1], fileName=fileName)
  
  if (nrow(filtsSig)>0){
    
    filtsSig$ID = filtsSigNames
    pcpDat <- melt(filtsSig[,c(1:7)], id.vars="ID")
    colnames(pcpDat) <- c("ID", "Sample", "Count")
    pcpDat$Sample <- as.character(pcpDat$Sample)
    
    plot_filteredSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[1], alpha=0.2) + xlab(paste("Cluster ", i, " (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + ylab("Count")
    
  } else{
    filtsSig = filtsSig[1,]
    filtsSig$ID = c("PlaceholderID")
    pcpDat <- melt(filtsSig[,c(1:7)], id.vars="ID")
    colnames(pcpDat) <- c("ID", "Sample", "Count")
    pcpDat$Sample <- as.character(pcpDat$Sample)
    plot_filteredSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[1], alpha=0) + xlab(paste("Cluster ", i, " (n=0", ")",sep="")) + ylab("Count")
  }
  
  jpeg(file = paste(getwd(), "/", outDir, "/", plotName, "_Sig_", nC, ".jpg", sep=""), width=1000, height=700)
  # We allow up to 4 plots in each column
  p = do.call("grid.arrange", c(append(plot_clustersSig, list(plot_filteredSig)), ncol=ceiling(nC/2)))
  invisible(dev.off())
}

for (i in 2:5){
  getPCP(i)
}
