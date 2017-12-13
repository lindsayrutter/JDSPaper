# This is the data from reference 6 of TMM Robinson paper

library(data.table)
load("LK_data.RData")
data = as.data.frame(MA.subsetA$M)
rownames(data) = as.character(MA.subsetA$genes$EnsemblGeneID)
setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","K.R1L1","L.R1L2","K.R1L3","L.R1L4","L.R1L6","K.R1L7","L.R1L8","K.R2L2","L.R2L3","K.R2L6")
data = as.data.frame(data)
data = data[,c(1,2,4,7,9,11,3,5,6,8,10)]

outDir = "DEG-tmm"

# Obtain R1 values
data <- data[,c(1:4,7:9)]

############# Get DEGs ############# 
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
plotDEG(rlogData2, metricList, outDir=outDir, option="scatterPrediction", threshVar="PValue", threshVal=1e-250)



# used for all except PI plots
#logData[,2:ncol(logData)] <- log(data[,2:ncol(logData)]+1)

plotDEG(logData, metricList, outDir=outDir, threshVar="PValue", threshVal=1e-250)
plotDEG(logData, metricList, outDir=outDir, option="volcano", threshVar="PValue", threshVal=1e-250)
plotDEG(logData, metricList, outDir=outDir, option="scatterOrthogonal", threshVar="PValue", threshVal=1e-250)
plotDEG(logData, metricList, outDir=outDir, option="scatterOrthogonal", threshVar="PValue", threshVal=1e-250, threshOrth=2)
plotDEG(logData, metricList, outDir=outDir, option="scatterOrthogonal", threshVar="PValue", threshVal=1e-250, threshOrth=1)
plotDEG(logData, metricList, outDir=outDir, option="scatterOrthogonal", threshVar="PValue", threshVal=1e-250, threshOrth=0.5)
plotDEG(logData, metricList, outDir=outDir, option="scatterPrediction", threshVar="PValue", threshVal=1e-250)
plotDEG(logData, metricList, outDir=outDir, option="parallelCoord", threshVar="PValue", threshVal=1e-250)

#NEW
for (nC in c(1:6)){plotClusters(data=logData, dataMetrics = metricList, nC=nC , threshVar = "PValue", threshVal = 1e-250, outDir=outDir)}

#Permutations
plotPermutations(data, nPerm = 10, topThresh = 100, outDir = outDir)
outDir = "DEG-tmm/log"
plotPermutations(logData, nPerm = 10, topThresh = 100, outDir = outDir)

# Statistics
length(which(lrt[[14]]$PValue<0.1/73320)) # 5704 had p-value <0.1/73320
length(which(lrt[[14]]$PValue<0.05/73320)) # 5546 had p-value <0.05/73320
length(which(lrt[[14]]$PValue<0.01/73320)) # 5224 had p-value <0.01/73320
length(which(lrt[[14]]$PValue<1e-250)) # 365 had p-value <0.01/73320


