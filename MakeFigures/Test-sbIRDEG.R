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

p1 = plotDEG(data, metricList, option="scatterPoints", threshVar="FDR", threshVal=0.05, pointSize=0.1)


dataMetrics=metricList; pointSize=0.5; degPointColor="orange"; threshVar="FDR"; threshVal=0.35



