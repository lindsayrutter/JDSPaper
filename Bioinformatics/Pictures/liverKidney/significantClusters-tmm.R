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

############# Get data metrics ############# 
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
sigL <- sigMets[which(sigMets$logFC<0),]
sigK <- sigMets[which(sigMets$logFC>0),]

RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

data_Rownames <- data$ID
data = data[,-1]
rownames(data) <- data_Rownames
#Normalize and log
cpm.data.new <- cpm(data, TRUE, TRUE)
# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(cpm.data.new, which="full", round=FALSE)
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

# File output information
plotName = "K_L"
outDir = "Clustering_data_FDR_05_TMM"

# Indices of the 9760 NAN rows. They had stdev=0 in the filt data
nID <- which(is.nan(dataqps$K.1))
# Set these filtered values that have all same values for samples to 0
dataqps[nID,1:6] <- 0

logSoy = origData
logSoy[,-1] <- log(origData[,-1]+1)

#####################################################

set.seed(1)
colList = scales::hue_pal()(4)
colList = colList[c(4,3,2,5,1)]

yMin = min(dataqps[,1:6])
yMax = max(dataqps[,1:6])

###########################

sigIDs = list(sigL$ID, sigK$ID)

plot_clustersSig = lapply(1:2, function(i){ 
  x = as.data.frame(dataqps[which(dataqps$ID %in% sigIDs[[i]]),])
  x$cluster = "color"
  x$cluster2 = factor(x$cluster)
  xNames = rownames(x)
  metricFDR = metrics[which(as.character(metrics$ID) %in% xNames),]
  sigID = metricFDR[metricFDR$FDR<0.05,]$ID
  xSig = x[which(rownames(x) %in% sigID),]
  xSigNames = rownames(xSig)
  nGenes = nrow(xSig)
  
  xSig$ID = xSigNames
  pcpDat <- melt(xSig[,c(1:7)], id.vars="ID")
  colnames(pcpDat) <- c("ID", "Sample", "Count")
  pcpDat$Sample <- as.character(pcpDat$Sample)
  
  pSig = ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[2], alpha=0.05) + ylab("Standardized Count") + ggtitle(paste("Significant Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=18), axis.text=element_text(size=18), axis.title=element_text(size=18))
  
  fileName = paste(getwd(), "/", outDir, "/", plotName, "_Sig_", i, ".jpg", sep="")
  jpeg(fileName)
  plot(pSig)
  invisible(dev.off())
  pSig
})

jpeg(file = paste(getwd(), "/", outDir, "/", plotName, "_Sig.jpg", sep=""), width=1000, height=700)
# We allow up to 4 plots in each column
p = do.call("grid.arrange", c(plot_clustersSig, ncol=2))
invisible(dev.off())

