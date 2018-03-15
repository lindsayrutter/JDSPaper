library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(GGally)
library(gridExtra)
library(cowplot)

data(soybean_cn)
data(soybean_cn_metrics)
data <- soybean_cn
metrics <- soybean_cn_metrics

# Focus on treatment groups S1 and S2
data <- data[,c(1:7)]

metricList = list()
metricList[["S1_S2"]] = metrics[["S1_S2"]]

p1 = plotDEG(data, metricList, option="scatterPoints", threshVar="FDR", threshVal=0.05, pointSize=0.5)
p2 = plotDEG(data, metricList, option="parallelCoord", threshVar="FDR", threshVal=0.05, lineSize = 0.5)
# default lineSize=0.1 in parallelCoord

plot1 <- plot_grid(ggmatrix_gtable(p1[["S1_S2"]]), labels=c("A"), ncol = 1, nrow = 1, label_size=9) + theme(plot.background = element_rect(size=0.1,linetype="solid",color="black"))
plot2 <- plot_grid(p2[["S1_S2"]], labels=c("B"), ncol = 1, nrow = 1, label_size=9) + theme(plot.background = element_rect(size=0.1,linetype="solid",color="black"))
grid.arrange(plot1, plot2, ncol=1, rel_heights = c(3/4, 1/4))
