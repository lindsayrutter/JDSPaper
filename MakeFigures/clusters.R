

source("functions.R")

data <- makeSbIRNoFiltObjects()[["soybean_ir_noFilt"]]
metrics <- makeSbIRNoFiltObjects()[["soybean_ir_noFilt_metrics"]][["N_P"]]
fulls <- makeSbIRNoFiltObjects()[["fulls"]]
dataqps <- makeSbIRNoFiltObjects()[["dataqps"]]

boxDat <- melt(fulls, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")

dendo = dataqps
rownames(dendo) = NULL
d = dist(as.matrix(dendo))
hc = hclust(d, method="ward.D")

# fileName = paste(getwd(), "/", outDir, "/", plotName, "_dendogram.jpg", sep="")
# jpeg(fileName)
# plot(hc, main="data Dendogram", xlab=NA, sub=NA)
# invisible(dev.off())

logSoy = soybean_ir
logSoy[,-1] <- log(soybean_ir[,-1]+1)






getPCP <- function(nC){
  
  set.seed(1)
  colList = scales::hue_pal()(nC+1)
  colList = colList[c(4,3,2,5,1)]
  k = cutree(hc, k=nC)
  
  yMin = min(dataqps[,1:6])
  yMax = max(dataqps[,1:6])
  
  ###########################
  
  # plot_clusters = lapply(1:nC, function(i){
  #   x = as.data.frame(dataqps[which(k==i),])
  #   nGenes = nrow(x)
  #   x$cluster = "color"
  #   x$cluster2 = factor(x$cluster)
  #   xNames = rownames(x)
  #   xFDR = metrics[which(metrics$ID %in% xNames),]$FDR
  #   
  #   x$ID = xNames
  #   
  #   pcpDat <- melt(x[,c(1:7)], id.vars="ID")
  #   colnames(pcpDat) <- c("ID", "Sample", "Count")
  #   boxDat$Sample <- as.character(boxDat$Sample)
  #   pcpDat$Sample <- as.character(pcpDat$Sample)
  #   
  #   p = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[i+1], alpha=0.05) + ylab("Standardized Count") + ggtitle(paste("Cluster ", i, " Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=18), axis.text=element_text(size=18), axis.title=element_text(size=18))
  #   
  #   fileName = paste(getwd(), "/", outDir, "/", plotName, "_", nC, "_", i, ".jpg", sep="")
  #   jpeg(fileName)
  #   plot(p)
  #   invisible(dev.off())
  #   p
  # })
  # ###########################
  # filts = as.data.frame(filts)
  # filts$cluster = "color"
  # filts$cluster2 = factor(filts$cluster)
  # nGenes = nrow(filts)
  # 
  # xNames = rownames(filts)
  # xFDR = metrics[which(metrics$ID %in% xNames),]$FDR
  # 
  # filts$ID = xNames
  # colnames(filts)[1:6] = colnames(dataqps)[1:6]
  # 
  # pcpDat <- melt(filts[,c(1:7)], id.vars="ID")
  # colnames(pcpDat) <- c("ID", "Sample", "Count")
  # pcpDat$Sample <- as.character(pcpDat$Sample)
  # 
  # plot_filtered = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[1], alpha=0.2) + ylab("Standardized Count") + ggtitle(paste("Filtered Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=18), axis.text=element_text(size=18), axis.title=element_text(size=18))
  # jpeg(file = paste(getwd(), "/", outDir, "/", plotName, "_", nC, ".jpg", sep=""), width=1000, height=700)
  # # We allow up to 4 plots in each column
  # p = do.call("grid.arrange", c(append(plot_clusters, list(plot_filtered)), ncol=ceiling(nC/2)))
  # invisible(dev.off())
  
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
    #saveRDS(xSigNames, file=paste0(getwd(), "/", outDir, "/Sig_", nC, "_", i, ".Rds"))
    
    if (nrow(xSig)>0){
      # scatMatMetrics = list()
      # scatMatMetrics[["N_P"]] = metrics[which(metrics$ID %in% xSigNames),]
      # scatMatMetrics[["N_P"]]$FDR = 10e-10
      # scatMatMetrics[["N_P"]]$ID = as.factor(as.character(scatMatMetrics[["N_P"]]$ID))
      # 
      # fileName = paste(getwd(), "/", outDir, "/", plotName, "_Sig_SM_", nC, "_", i, ".jpg", sep="")
      # ret <- plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "FDR", threshVal = 0.05, degPointColor = colList[i+1], fileName=fileName)
      # jpeg(fileName, height=700, width=700)
      # ret[[plotName]] + xlab("Logged Count") + ylab("Logged Count") + ggtitle(paste("Cluster ", i, " Significant Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=14), axis.text=element_text(size=14), axis.title=element_text(size=18), strip.text = element_text(size = 14))
      # invisible(dev.off())
      
      xSig$ID = xSigNames
      pcpDat <- melt(xSig[,c(1:7)], id.vars="ID")
      colnames(pcpDat) <- c("ID", "Sample", "Count")
      pcpDat$Sample <- as.character(pcpDat$Sample)
      
      pSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[i+1], alpha=0.5) + ylab("Standardized Count") + ggtitle(paste("Cluster ", i, " Significant Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=18), axis.text=element_text(size=18), axis.title=element_text(size=18))
    }
    #else{
    #   plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "FDR", threshVal = 0.05, degPointColor = colList[i+1], fileName=fileName)
    #   ret <- plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "FDR", threshVal = 0.05, degPointColor = colList[i+1], fileName=fileName)
    #   fileName = paste(getwd(), "/", outDir, "/", plotName, "_Sig_SM_", nC, "_", i, ".jpg", sep="")
    #   jpeg(fileName)
    #   ret[[plotName]] + xlab("Logged Count") + ylab("Logged Count") + ggtitle(paste("Cluster ", i, " Significant Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5))
    #   invisible(dev.off())
    #   
    #   pSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + ylab("Standardized Count") + ggtitle(paste("Cluster ", i, " Significant Genes (n=0)")) + theme(plot.title = element_text(hjust = 0.5, size=14), axis.text=element_text(size=14), axis.title=element_text(size=14))
    # }
    #fileName = paste(getwd(), "/", outDir, "/", plotName, "_Sig_", nC, "_", i, ".jpg", sep="")
    #jpeg(fileName)
    #plot(pSig)
    #invisible(dev.off())
    pSig
  })
  
  jpeg(file = paste(getwd(), "/", outDir, "/", plotName, "_Sig_", nC, ".jpg", sep=""), width=1000, height=700)
  # We allow up to 4 plots in each column
  p = do.call("grid.arrange", c(plot_clustersSig, ncol=ceiling(nC/2)))
  invisible(dev.off())
}