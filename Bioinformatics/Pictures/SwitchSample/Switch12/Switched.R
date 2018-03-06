library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)

data(soybean_cn)
load("../../soybean_cn_switch12_metrics.rda")
data <- soybean_cn
metrics <- soybean_cn_switch12_metrics

# Only focus on S1_S2
data <- data[,c(1:7)]
data <- data[,c(1,2,3,5,4,6,7)]
colnames(data) <- c("ID","S1.1","S1.2","S1.3","S2.1","S2.2","S2.3")
baseOutDir = "/Users/lindz/JDSPaper/Bioinformatics/Pictures/SwitchSample/Switch12"

# Obtain Group values
outDir = paste0(baseOutDir, "/Switched")
dataSel <- data
plotScatterStatic(dataSel, outDir = outDir, option ="point")
plotScatterStatic(dataSel, outDir = outDir)
boxSel = dataSel[,-1] %>% gather(Sample,Count)
boxSel$group = c(rep("S1",nrow(boxSel)/2), rep("S2", nrow(boxSel)/2))
bPlot = ggplot(boxSel, aes(x=Sample, y=Count, fill=group)) + geom_boxplot() + scale_fill_manual(values=c("limegreen","magenta")) + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")
png(filename=paste0(outDir = outDir,"/boxplot.jpg"), width=240, height=240)
bPlot
dev.off()
plotScatterStatic(dataSel, threshOrth = 4, outDir = outDir, option="orthogonal")
plotScatterStatic(dataSel, piLevel=0.99999, outDir = outDir, option="prediction")

colors <- c(rep(c("limegreen"), 3), rep(c("magenta"), 3))
#p <- plotMDS(dataSel[,2:ncol(dataSel)], col=colors, xlab="Dim 1", ylab="Dim 2", labels=colnames(dataSel[,2:ncol(dataSel)]))

plotMDS(dataSel[,2:ncol(dataSel)], col=colors, xlab="Dim 1", ylab="Dim 2", labels=colnames(dataSel[,2:ncol(dataSel)]), xlim=c(-4, 4), ylim=c(-2, 2.5))

tDat <- t(dataSel[,2:ncol(dataSel)])
datD <- as.matrix(dist(tDat))
fit <- cmdscale(datD, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
qplot(x,y) + geom_text(label = names(x), nudge_y = 0.35, color = 'royalblue', fontface="bold") + labs(x = "Dim 1", y = "Dim 2") + theme(text = element_text(size=12))



plotScatterInteractive(dataSel, xbins=20)

############# Get DEGs ############# 
metricList = list()
metricList[["S1_S2"]] = metrics[["S1_S2"]]

############# Create DEG plots #############
plotDEG(data, metricList, outDir=outDir, threshVar="FDR", threshVal=0.05)
plotDEG(data, metricList, outDir=outDir, option="volcano", threshVar="FDR", threshVal=0.05)
plotDEG(data, metricList, outDir=outDir, option="scatterOrthogonal", threshVar="FDR", threshVal=0.05)
plotDEG(data, metricList, outDir=outDir, option="scatterPrediction", threshVar="FDR", threshVal=0.05)
plotDEG(data, metricList, outDir=outDir, option="parallelCoord", threshVar="FDR", threshVal=0.05)

plotRepLine(data, metricList)
