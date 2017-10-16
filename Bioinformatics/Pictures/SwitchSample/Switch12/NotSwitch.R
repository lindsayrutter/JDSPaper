library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)

data(soybean_cn)
data(soybean_cn_metrics)
data <- soybean_cn
metrics <- soybean_cn_metrics

# Only focus on S1_S2
data <- data[,c(1:7)]
baseOutDir = "/Users/lindz/JDSPaper/Bioinformatics/Pictures/SwitchSample/Switch12"

# Obtain Group values
outDir = paste0(baseOutDir, "/NotSwitch")
dataSel <- data
plotScatterStatic(dataSel, outDir = outDir, option="point")
plotScatterStatic(dataSel, outDir = outDir, option="point", pointSize=0.5)
plotScatterStatic(dataSel, outDir = outDir, option="point", pointSize=0.25)
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
metricList = list()
metricList[["S1_S2"]] = metrics[["S1_S2"]]

############# Create DEG plots #############
plotDEG(data, metricList, outDir=outDir, threshVar="FDR", threshVal=0.05)
plotDEG(data, metricList, option="scatterPoints", outDir=outDir, threshVar="FDR", threshVal=0.05, pointSize=0.5)
plotDEG(data, metricList, outDir=outDir, option="volcano", threshVar="FDR", threshVal=0.05)
plotDEG(data, metricList, outDir=outDir, option="scatterOrthogonal", threshVar="FDR", threshVal=0.05)
plotDEG(data, metricList, outDir=outDir, option="scatterPrediction", threshVar="FDR", threshVal=0.05)
plotDEG(data, metricList, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.05)

plotRepLine(data, metricList)

