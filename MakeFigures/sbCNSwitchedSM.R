library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)

data(soybean_cn)
load("sbCNSwitchedMetrics.rda")
data <- soybean_cn
metrics <- sbCNSwitchedMetrics

# Focus on treatment groups S1 and S2
data <- data[,c(1:7)]
data <- data[,c(1,2,3,5,4,6,7)]
colnames(data) <- c("ID","S1.1","S1.2","S1.3","S2.1","S2.2","S2.3")

p = plotScatterStatic(data, option="point", pointSize=0.25, saveFile = FALSE)

# Plot figure
p[["S1_S2"]]
