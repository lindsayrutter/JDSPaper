library(bigPint)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(gridExtra)

data(soybean_ir)
data <- soybean_ir

# Focus on treatment groups S1 and S2
data <- data[,c(1:7)]
# MDS plot of non-switched data
colors <- c(rep(c("blueviolet"), 3), rep(c("chocolate1"), 3))
ggplot(plotMDS(data[,2:ncol(data)], col=colors, xlab="Dim 1", ylab="Dim 2", labels=colnames(data[,2:ncol(data)]), xlim=c(-2800, 2800), ylim=c(-2800, 2800)))

# Logged
data[,c(2:7)] <- log(data[,c(2:7)])
ggplot(plotMDS(data[,2:ncol(data)], col=colors, xlab="Dim 1", ylab="Dim 2", labels=colnames(data[,2:ncol(data)]), xlim=c(-1.6, 1.6), ylim=c(-1.6, 1.6)))

