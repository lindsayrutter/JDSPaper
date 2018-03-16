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
p2 = plotDEG(data, metricList, option="parallelCoord", threshVar="FDR", threshVal=0.05, lineSize = 0.3)

plot1 <- ggmatrix_gtable(p1[["S1_S2"]])
plot2 <- p2[["S1_S2"]] + theme(axis.text=element_text(size=10), axis.title =element_text(size=10))
plot_grid(plot1, plot2, align = "v", nrow = 2, rel_heights = c(3/4, 1/4), labels = c("A", "B"), label_size = 9)

