library(edgeR)
library(DESeq2)
library(bigPint)
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
rlogData <- rlog(as.matrix(logData[,2:7]))
rlogData <- as.data.frame(rlogData)
rlogData2 <- cbind(ID = data$ID, rlogData)
rlogData2$ID <- as.character(rlogData2$ID)

# There are 40 points that intersect all blue points and union all red points from all 15 subplots
redBluePoints <- plotDEG(data=rlogData2, dataMetrics=metricList, outDir=outDir, option="scatterPrediction", threshVar="PValue", threshVal=1e-250, piLevel=0.95)
redBluePoints[["K_L"]] <- as.character(redBluePoints[["K_L"]]$ID[which(redBluePoints[["K_L"]]$Color=="Red")])
plotDEG(rlogData2, lineList = redBluePoints, lineSize=0.5, lineColor = "red", outDir=outDir, option ="parallelCoord") # Change name to K_L_deg_pcp_0.05_RED.jpg

redBluePoints <- plotDEG(data=rlogData2, dataMetrics=metricList, outDir=outDir, option="scatterPrediction", threshVar="PValue", threshVal=1e-250, piLevel=0.95)
redBluePoints[["K_L"]] <- as.character(redBluePoints[["K_L"]]$ID[which(redBluePoints[["K_L"]]$Color=="Blue")])
plotDEG(rlogData2, lineList = redBluePoints, lineSize=0.5, lineColor = "blue", outDir=outDir, option ="parallelCoord") # Change name to K_L_deg_pcp_0.05_BLUE.jpg

