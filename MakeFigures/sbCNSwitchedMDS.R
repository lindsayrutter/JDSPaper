library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(gridExtra)

data(soybean_cn)
load("sbCNSwitchedMetrics.rda")
data <- soybean_cn
metrics <- sbCNSwitchedMetrics

# Focus on treatment groups S1 and S2
data <- data[,c(1:7)]
dataSwitched <- data[,c(1,2,3,5,4,6,7)]
colnames(dataSwitched) <- c("ID","S1.1","S1.2","S1.3","S2.1","S2.2","S2.3")

# Boxplot of non-switched data
boxSel = data[,-1] %>% gather(Sample,Count)
boxSel$group = c(rep("S1",nrow(boxSel)/2), rep("S2", nrow(boxSel)/2))
bPlot = ggplot(boxSel, aes(x=Sample, y=Count, fill=group)) + geom_boxplot() + scale_fill_manual(values=c("limegreen","magenta")) + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

# Boxplot of switched data
boxSel = dataSwitched[,-1] %>% gather(Sample,Count)
boxSel$group = c(rep("S1",nrow(boxSel)/2), rep("S2", nrow(boxSel)/2))
bPlotSwitched = ggplot(boxSel, aes(x=Sample, y=Count, fill=group)) + geom_boxplot() + scale_fill_manual(values=c("limegreen","magenta")) + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

# MDS plot of non-switched data
colors <- c(rep(c("limegreen"), 3), rep(c("magenta"), 3))
ggplot(plotMDS(data[,2:ncol(data)], col=colors, xlab="Dim 1", ylab="Dim 2", labels=colnames(data[,2:ncol(data)]), xlim=c(-4, 4), ylim=c(-4, 4)))
ggplot(plotMDS(data[,2:ncol(data)], col=colors, xlab="Dim 1", ylab="Dim 2", pch=20, xlim=c(-4, 4), ylim=c(-4, 4)))

mdsDat = data[,2:ncol(data)]
ggplot(mdsDat, aes(x,y)) + geom_text(data = mdsDat[c(1:3),], label = rownames(mdsDat[c(1:3),]), nudge_y = 0.35, fontface="bold", color = "royalblue") + geom_text(data = mdsDat[c(4:6),], label = rownames(mdsDat[c(4:6),]), nudge_y = 0.35, fontface="bold", color = "darkorange2") + labs(x = "Dim 1", y = "Dim 2") + theme(text = element_text(size=12))



# MDS plot of switched data
plotMDS(dataSwitched[,2:ncol(dataSwitched)], col=colors, xlab="Dim 1", ylab="Dim 2", labels=colnames(dataSwitched[,2:ncol(dataSwitched)]), xlim=c(-4, 4), ylim=c(-4, 4))
plotMDS(dataSwitched[,2:ncol(dataSwitched)], col=colors, xlab="Dim 1", ylab="Dim 2", pch=20, xlim=c(-4, 4), ylim=c(-4, 4))



# Plot all clusters as a grid
do.call("grid.arrange", c(append(plot_clusters, list(plot_filtered)), ncol=ceiling(nC/2)))

grid.arrange(bPlot, bPlotSwitched, plotMDS(data[,2:ncol(data)], col=colors, xlab="Dim 1", ylab="Dim 2", labels=colnames(data[,2:ncol(data)]), xlim=c(-4, 4), ylim=c(-2, 2.5))
, plotMDS(dataSwitched[,2:ncol(dataSwitched)], col=colors, xlab="Dim 1", ylab="Dim 2", labels=colnames(dataSwitched[,2:ncol(dataSwitched)]), xlim=c(-4, 4), ylim=c(-2, 2.5))
, ncol=2)
