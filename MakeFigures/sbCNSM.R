library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)

data(soybean_cn)
data(soybean_cn_metrics)
data <- soybean_cn
metrics <- soybean_cn_metrics

# Focus on treatment groups S1 and S2
data <- data[,c(1:7)]
p = plotScatterStatic(data, option="point", pointSize=0.25, saveFile = FALSE)

# Plot figure
p[["S1_S2"]]
