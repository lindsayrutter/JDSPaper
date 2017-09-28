library(bigPint)

data(soybean_cn)
data <- soybean_cn
data <- data[,c(1:7)]
baseOutDir = "/Users/lindz/JDSPaper/SwitchSample"

# Obtain Group values
outDir = paste0(baseOutDir, "/NotSwitch")
dataSel <- data
plotScatterStatic(dataSel, outDir = outDir)
boxSel = dataSel[,-1] %>% gather(Sample,Count)
bPlot = ggplot(boxSel, aes(x=Sample, y=Count)) + geom_boxplot()
png(filename=paste0(outDir = outDir,"/boxplot.jpg"))
bPlot
dev.off()
plotScatterStatic(dataSel, threshOrth = 4, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, piLevel=0.99999, outDir = outDir, option="prediction")
plotMDS(dataSel[,2:ncol(dataSel)])
plotScatterInteractive(dataSel, xbins=20)

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
metricList[["S1_S2"]] = ret

############# Create DEG plots #############
plotDEG(data, metricList, outDir=outDir, threshVar="PValue", threshVal=0.001)
plotDEG(data, metricList, outDir=outDir, option="volcano", threshVar="PValue", threshVal=0.001)
plotDEG(data, metricList, outDir=outDir, option="scatterOrthogonal", threshVar="PValue", threshVal=0.001)
plotDEG(data, metricList, outDir=outDir, option="scatterPrediction", threshVar="PValue", threshVal=0.001)
plotDEG(data, metricList, outDir=outDir, option="parallelCoord", threshVar="PValue", threshVal=0.001)

plotRepLine(data, metricList)

# Statistics
length(which(lrt[[14]]$PValue<0.1/73320)) # 13 had p-value <0.1/73320
length(which(lrt[[14]]$PValue<0.05/73320)) # 7 had p-value <0.05/73320
length(which(lrt[[14]]$PValue<0.01/73320)) # 5 had p-value <0.01/73320
length(which(lrt[[14]]$PValue<0.001)) # 368 had p-value <0.001
length(which(lrt[[14]]$PValue<0.01)) # 1164 had p-value <0.01