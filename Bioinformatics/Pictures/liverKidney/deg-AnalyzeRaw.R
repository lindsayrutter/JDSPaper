# This is the data from reference 6 of TMM Robinson paper

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
ret = data.frame(ID=rownames(dataSel), lrt[[14]])
ret$ID = as.character(ret$ID)
ret = as.data.frame(ret)
metricList = list()
metricList[["K_L"]] = ret

############# Create DEG plots #############
logData <- data
logData[,2:ncol(logData)] <- log(data[,2:ncol(logData)]+1)
plotDEG(logData, metricList, outDir=outDir, option="scatterPrediction", threshVar="PValue", threshVal=1e-250)

# Statistics
length(which(lrt[[14]]$PValue<0.1/73320)) # 7040 had p-value <0.1/73320
length(which(lrt[[14]]$PValue<0.05/73320)) # 6870 had p-value <0.05/73320
length(which(lrt[[14]]$PValue<0.01/73320)) # 6491 had p-value <0.01/73320
length(which(lrt[[14]]$PValue<1e-250)) # 358 had p-value <1e-250 (most of these are pretty much 0 and cannot be ordered)



