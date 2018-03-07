# This is the data from reference 6 of TMM Robinson paper

library(edgeR)
library(DESeq2)
library(bigPint)
library(data.table)

load("LK_data.RData")
data = as.data.frame(MA.subsetA$M)
rownames(data) = as.character(MA.subsetA$genes$EnsemblGeneID)
setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","K.R1L1","L.R1L2","K.R1L3","L.R1L4","L.R1L6","K.R1L7","L.R1L8","K.R2L2","L.R2L3","K.R2L6")
data = as.data.frame(data)
data = data[,c(1,2,4,7,9,11,3,5,6,8,10)]

outDir = "DEG-raw"

# Obtain R1 values
data <- data[,c(1:4,7:9)]
colnames(data) <- c("ID", "K.1", "K.2", "K.3", "L.1", "L.2", "L.3")

############# Get DEGs ############# 
rowNames = data[,1]
dataSel = data[,-1]
rownames(dataSel) = rowNames
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=dataSel,group=group)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
ret = data.frame(ID=rownames(topTags(lrt, n = nrow(y[[1]]))[[1]]), topTags(lrt, n = nrow(y[[1]]))[[1]])
ret$ID = as.character(ret$ID)
ret = as.data.frame(ret)
metricList = list()
metricList[["K_L"]] = ret

metricList2 = metricList
metricList2[["K_L"]]$FDR = 1


############# Create DEG plots #############
logData <- data
rownames(logData) <- logData$ID
logData[,2:ncol(logData)] <- log(data[,2:ncol(logData)]+1)
#plotDEG(logData, metricList, outDir=outDir, option="scatterPrediction", threshVar="FDR", threshVal=0.0001)

ret <- plotDEG(data = logData, dataMetrics = metricList2, outDir=outDir, option="scatterPoints", threshVar = "FDR", threshVal=1e-3)
retP <- ret[["K_L"]]
fileName = paste(getwd(), "/", outDir, "/K_L_SM.jpg", sep="")
jpeg(fileName, height=700, width=700)
retP + xlab("Logged Count") + ylab("Logged Count") +theme(plot.title = element_text(hjust = 0.5, size=15), axis.text=element_text(size=14), axis.title=element_text(size=15), strip.text = element_text(size = 13))
dev.off()

ret <- plotDEG(data = logData, dataMetrics = metricList, outDir=outDir, option="scatterPrediction", threshVar = "FDR", threshVal=1e-12)
retP <- ret[["K_L"]]
fileName = paste(getwd(), "/", outDir, "/K_L_PI_1e-12.jpg", sep="")
jpeg(fileName, height=700, width=700)
retP + xlab("Logged Count") + ylab("Logged Count") + ggtitle(paste("Significant Genes (n=", format(length(which(metricList[["K_L"]]$FDR<1e-12)), big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=15), axis.text=element_text(size=14), axis.title=element_text(size=15), strip.text = element_text(size = 13))
dev.off()

ret <- plotDEG(data = logData, dataMetrics = metricList, outDir=outDir, option="scatterPrediction", threshVar = "FDR", threshVal=1e-21)
retP <- ret[["K_L"]]
fileName = paste(getwd(), "/", outDir, "/K_L_PI_1e-21.jpg", sep="")
jpeg(fileName, height=700, width=700)
retP + xlab("Logged Count") + ylab("Logged Count") + ggtitle(paste("Significant Genes (n=", format(length(which(metricList[["K_L"]]$FDR<1e-21)), big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=15), axis.text=element_text(size=14), axis.title=element_text(size=15), strip.text = element_text(size = 13))
dev.off()

ret <- plotDEG(data = logData, dataMetrics = metricList, outDir=outDir, option="scatterPrediction", threshVar = "FDR", threshVal=1e-39)
retP <- ret[["K_L"]]
fileName = paste(getwd(), "/", outDir, "/K_L_PI_1e-39.jpg", sep="")
jpeg(fileName, height=700, width=700)
retP + xlab("Logged Count") + ylab("Logged Count") + ggtitle(paste("Significant Genes (n=", format(length(which(metricList[["K_L"]]$FDR<1e-39)), big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=15), axis.text=element_text(size=14), axis.title=element_text(size=15), strip.text = element_text(size = 13))
dev.off()

ret <- plotDEG(data = logData, dataMetrics = metricList, outDir=outDir, option="scatterPoints", threshVar = "FDR", threshVal=1e-39, pointSize=0.001)
retP <- ret[["K_L"]]
fileName = paste(getwd(), "/", outDir, "/K_L_SM_1e-39_PointSize001.jpg", sep="")
jpeg(fileName, height=700, width=700)
retP + xlab("Logged Count") + ylab("Logged Count") + ggtitle(paste("Significant Genes (n=", format(length(which(metricList[["K_L"]]$FDR<1e-39)), big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=15), axis.text=element_text(size=14), axis.title=element_text(size=15), strip.text = element_text(size = 13))
dev.off()

# Statistics
length(which(ret$FDR<0.001))
