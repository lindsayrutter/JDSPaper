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
library(data.table)

############# Get data ############# 
load("LK_data.RData")
data = as.data.frame(MA.subsetA$M)
rownames(data) = as.character(MA.subsetA$genes$EnsemblGeneID)
setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","K.R1L1","L.R1L2","K.R1L3","L.R1L4","L.R1L6","K.R1L7","L.R1L8","K.R2L2","L.R2L3","K.R2L6")
data = as.data.frame(data)
data = data[,c(1,2,4,7,9,11,3,5,6,8,10)]
# Obtain R1 values
data <- data[,c(1:4,7:9)]
colnames(data) <- c("ID", "K.1", "K.2", "K.3", "L.1", "L.2", "L.3")
origData = data

############# Get raw data metrics ############# 
rowNames = data[,1]
dataSel = data[,-1]
rownames(dataSel) = rowNames
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=dataSel,group=group)
#y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
ret = data.frame(ID=rownames(topTags(lrt, n = nrow(y[[1]]))[[1]]), topTags(lrt, n = nrow(y[[1]]))[[1]])
ret$ID = as.character(ret$ID)
ret = as.data.frame(ret)
metricList = list()
metricList[["K_L"]] = ret
metrics <- metricList[["K_L"]]
sigMets = metrics[which(metrics$FDR<0.001),]
sigK_Raw <- sigMets[which(sigMets$logFC<0),]
sigL_Raw <- sigMets[which(sigMets$logFC>0),]

############# Get TMM data metrics ############# 
rowNames = data[,1]
dataSel = data[,-1]
rownames(dataSel) = rowNames
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=dataSel,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
ret = data.frame(ID=rownames(topTags(lrt, n = nrow(y[[1]]))[[1]]), topTags(lrt, n = nrow(y[[1]]))[[1]])
ret$ID = as.character(ret$ID)
ret = as.data.frame(ret)
metricList = list()
metricList[["K_L"]] = ret
metrics <- metricList[["K_L"]]
sigMets = metrics[which(metrics$FDR<0.001),]
sigK_TMM <- sigMets[which(sigMets$logFC<0),]
sigL_TMM <- sigMets[which(sigMets$logFC>0),]

length(which(sigL_TMM$ID %in% sigL_Raw$ID)) # All 1968 genes in Raw liver were in TMM liver
length(which(sigL_Raw$ID %in% sigL_TMM$ID)) # 1968 genes in Raw liver remain in TMM liver
length(which(!sigL_TMM$ID %in% sigL_Raw$ID)) # 1578 genes that were not in Raw liver were introduced in TMM liver

origDEG <- sigL_Raw


# File output information
plotName = "K_L"
outDir = "Clustering_data_FDR_001_TMMvRaw_Orig"
logSoy = origData
logSoy[,-1] <- log(origData[,-1]+1)
totalColor = scales::seq_gradient_pal("darkorange", "orangered3", "Lab")(seq(0,1,length.out=8))[4]

# Make total scatterplot matrix
scatMatMetrics = list()
scatMatMetrics[["K_L"]] = metrics[which(metrics$ID %in% origDEG$ID),]
scatMatMetrics[["K_L"]]$FDR = 10e-10
scatMatMetrics[["K_L"]]$ID = as.factor(as.character(scatMatMetrics[["K_L"]]$ID))
fileName = paste(getwd(), "/", outDir, "/", plotName, "_Sig_SM_Orig.jpg", sep="")
ret <- plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "FDR", threshVal = 0.05, degPointColor = totalColor, fileName=fileName)
jpeg(fileName, height=700, width=700)
ret[[plotName]] + xlab("Logged Count") + ylab("Logged Count") + ggtitle(paste("Significant Genes Original (n=", format(nrow(scatMatMetrics[["K_L"]]), big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=14), axis.text=element_text(size=14), axis.title=element_text(size=18), strip.text = element_text(size = 14))
invisible(dev.off())


RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

data_Rownames <- data$ID
data = data[,-1]
rownames(data) <- data_Rownames
#Normalize and log (This allow standardized medians to be closer to 0)
#cpm.data.new <- cpm(data, TRUE, TRUE)
# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(as.matrix(data), which="full", round=FALSE)
data = as.data.frame(data)
# Add mean and standard deviation for each row/gene
data = mutate(data, mean = (K.1+K.2+K.3+L.1+L.2+L.3)/6, stdev = RowSD(cbind(K.1,K.2,K.3,L.1,L.2,L.3)))
rownames(data)=data_Rownames
data$ID <- data_Rownames
# Remove the genes with lowest quartile of mean and standard deviation

dataqps <- t(apply(as.matrix(data[,1:6]), 1, scale))
dataqps <- as.data.frame(dataqps)
colnames(dataqps) <- colnames(data[,1:6])
dataqps$ID <- rownames(dataqps)

# Comine the filtered and remaining data
fulls <- dataqps
boxDat <- melt(fulls, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")

# Indices of the 775 NAN rows. They had stdev=0 in the filt data
nID <- which(is.nan(dataqps$K.1))
# Set these filtered values that have all same values for samples to 0
dataqps[nID,1:6] <- 0

#####################################################

yMin = min(dataqps[,1:6])
yMax = max(dataqps[,1:6])

###########################

x = as.data.frame(dataqps[which(dataqps$ID %in% origDEG$ID),])
x$cluster = "color"
x$cluster2 = factor(x$cluster)
xNames = rownames(x)
xSig = x
xSigNames = rownames(xSig)
nGenes = nrow(xSig)

xSig$ID = xSigNames
pcpDat <- melt(xSig[,c(1:7)], id.vars="ID")
colnames(pcpDat) <- c("ID", "Sample", "Count")
pcpDat$Sample <- as.character(pcpDat$Sample)

pSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = "darkmagenta", alpha=0.1) + ylab("Standardized Count") + theme(plot.title = element_text(hjust = 0.5, size=25), axis.text=element_text(size=25), axis.title=element_text(size=25)) + ggtitle(paste("Significant Liver Genes Raw Original (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep=""))

fileName = paste(getwd(), "/", outDir, "/", plotName, "_KSigOrigWithTMM.jpg", sep="")
jpeg(fileName)
plot(pSig)
invisible(dev.off())



######## Cluster 3974 significant ones ##########
dendo = xSig[,1:6] # or dataqps? (If do fulls, then have NAs introduced by conversion)
rownames(dendo) = NULL
d = dist(as.matrix(dendo))
hc = hclust(d, method="ward.D")

fileName = paste(getwd(), "/", outDir, "/", plotName, "_dendogram.jpg", sep="")
jpeg(fileName)
plot(hc, main="data Dendogram", xlab=NA, sub=NA)
invisible(dev.off())


getPCP <- function(nC){
  
  set.seed(1)
  colList = scales::seq_gradient_pal("darkorange", "orangered3", "Lab")(seq(0,1,length.out=nC))
  k = cutree(hc, k=nC)
  
  plot_clusters = lapply(1:nC, function(j){
    i = rev(order(table(k)))[j]
    x = as.data.frame(xSig[,1:6][which(k==i),])
    nGenes = nrow(x)
    x$cluster = "color"
    x$cluster2 = factor(x$cluster)
    xNames = rownames(x)
    x$ID = xNames
    xSigNames = rownames(x)
    saveRDS(xSigNames, file=paste0(getwd(), "/", outDir, "/Sig_", nC, "_", j, ".Rds"))
    
    
    
    scatMatMetrics = list()
    scatMatMetrics[["K_L"]] = metrics[which(metrics$ID %in% x$ID),]
    scatMatMetrics[["K_L"]]$FDR = 10e-10
    scatMatMetrics[["K_L"]]$ID = as.factor(as.character(scatMatMetrics[["K_L"]]$ID))
    
    fileName = paste(getwd(), "/", outDir, "/", plotName, "_Sig_SM_Orig_", nC, "_", j, ".jpg", sep="")
    ret <- plotDEG(data = logSoy, dataMetrics = scatMatMetrics, option="scatterPoints", threshVar = "FDR", threshVal = 0.05, degPointColor = colList[j], fileName=fileName)
    jpeg(fileName, height=700, width=700)
    ret[[plotName]] + xlab("Logged Count") + ylab("Logged Count") + ggtitle(paste("Cluster ", j, " Significant Original Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=14), axis.text=element_text(size=14), axis.title=element_text(size=18), strip.text = element_text(size = 14))
    invisible(dev.off())
    
    
    
    pcpDat <- melt(x[,c(1:6,9)], id.vars="ID")
    colnames(pcpDat) <- c("ID", "Sample", "Count")
    boxDat$Sample <- as.character(boxDat$Sample)
    pcpDat$Sample <- as.character(pcpDat$Sample)
    
    p = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[j], alpha=0.5) + ylab("Standardized Count") + ggtitle(paste("Cluster ", j, " Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=32), axis.text=element_text(size=32), axis.title=element_text(size=32))
    
    fileName = paste(getwd(), "/", outDir, "/", plotName, "_", nC, "_", i, ".jpg", sep="")
    jpeg(fileName, width=650, height=500)
    plot(p)
    invisible(dev.off())
    p
  })
  jpeg(file = paste(getwd(), "/", outDir, "/", plotName, "_", nC, ".jpg", sep=""), width=1000, height=1800)
  # We allow up to 4 plots in each column
  if (nC %in% c(2,3,4)){
    p = do.call("grid.arrange", c(plot_clusters, ncol=2))
  }
  else{
    p = do.call("grid.arrange", c(plot_clusters, ncol=ceiling(nC/4)))
  }
  invisible(dev.off())
}

for (i in 2:8){
  getPCP(i)
}
