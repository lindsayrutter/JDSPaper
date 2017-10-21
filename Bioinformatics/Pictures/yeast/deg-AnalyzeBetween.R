# This is the data from EDASeq vignette

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(edgeR)
data = readRDS("dataBetween.rds")
setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","Y1.1","Y1.2","Y2.1","Y2.2","Y7.1","Y7.2","Y4.1","Y4.2","D.1","D.2","D.7","G.1","G.2","G.3")
data = as.data.frame(data)

outDir = "/Users/lindz/JDSPaper/Data/risso-Yeast/DEG-between"

############# Get DEGs ############# 
rowNames = data[,1]
dataSel = data[,-1]
rownames(dataSel) = rowNames
group <- factor(c("Y1","Y1","Y2","Y2","Y7","Y7","Y4","Y4","D","D","D","G","G","G"))
y <- DGEList(counts=dataSel,group=group)
design <- model.matrix(~group)
colnames(design) = levels(group)
y <- estimateDisp(y, design)
fit <- glmFit(y,design)

metricList <- list()
for (i in 1:(ncol(fit)-1)){
  for (j in (i+1):ncol(fit)){
    contrast=rep(0,ncol(fit))
    contrast[i]=1
    contrast[j]=-1
    lrt <- glmLRT(fit, contrast=contrast)
    lrt <- topTags(lrt, n = nrow(y[[1]]))[[1]]
    
    setDT(lrt, keep.rownames = TRUE)[]
    colnames(lrt)[1] = "ID"
    lrt <- as.data.frame(lrt)
    
    metricList[[paste0(colnames(fit)[i], "_", colnames(fit)[j])]] <- lrt
  }
}
betweenMetricList = metricList
saveRDS(betweenMetricList, "metricListBetween.rds")

############# Create DEG plots #############
plotData <- data
plotData[,2:ncol(plotData)] <- log(data[,2:ncol(plotData)]+1)

plotDEG(plotData, metricList, outDir=outDir, threshVar="FDR", threshVal=0.05)
plotDEG(plotData, metricList, outDir=outDir, threshVar="FDR", threshVal=0.01)

plotDEG(plotData, metricList, outDir=outDir, option = "volcano", threshVar="FDR", threshVal=0.05)
plotDEG(plotData, metricList, outDir=outDir, option = "volcano", threshVar="FDR", threshVal=0.01)

plotDEG(plotData, metricList, outDir=outDir, option="scatterOrthogonal", threshVar="FDR", threshOrth = 1, threshVal=0.05)
plotDEG(plotData, metricList, outDir=outDir, option="scatterOrthogonal", threshVar="FDR", threshOrth = 1, threshVal=0.01)

plotDEG(plotData, metricList, outDir=outDir, option="scatterPrediction", threshVar="FDR", threshVal=0.05)
plotDEG(plotData, metricList, outDir=outDir, option="scatterPrediction", threshVar="FDR", threshVal=0.01)

plotDEG(plotData, metricList, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.05, lineSize = 0.3)
plotDEG(plotData, metricList, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.01, lineSize = 0.3)

plotRepLine(plotData, metricList)
