library(rtracklayer)
library(Rsamtools)
library(grid)
library(GenomicAlignments)
library(ggplot2)
library(GGally)
library(edgeR)
library(stringr)
library(EDASeq)
library(dplyr)
library(matrixStats)
library(gridExtra)
library(reshape2)
library(scales)
library(bigPint)

source("functions.R")

data("soybean_ir")
#data("soybean_ir_metrics")
data <- soybean_ir
#metrics <- soybean_ir_metrics[["N_P"]]
load("../Bioinformatics/Pictures/FilterNotSig/soybean_ir_noFilt_metrics.rda")
metrics <- soybean_ir_noFilt_metrics[["N_P"]]

# Filter, normalize, and standardize the data so each gene has mean=0 and stdev=1
res <- filterStandardize(data)
# Fitered data standardized
filts <- res[["filts"]]
# Non-filtered data standardized
datas <- res[["datas"]]
# Hierarchical clustering object
hc <- res[["hc"]]
# Full data standardized
fulls <- rbind(datas, filts)

# Number of clusters
nC = 4
# Number of samples
nCol = 6
colList = scales::hue_pal()(nC+1)
colList <- colList[c(4, 3, 2, 5, 1)]
# Hierarchical clustering
k = cutree(hc, k=nC)

yMin = min(datas[,1:nCol])
yMax = max(datas[,1:nCol])

# Create background boxplot data
boxDat <- melt(fulls, id.vars="ID")
colnames(boxDat) <- c("ID", "Sample", "Count")

logSoy = soybean_ir
logSoy[,-1] <- log(soybean_ir[,-1]+1)

# Plot scatterplot matrix for Cluster 1 significant genes
plotClusterSM(1)
# Plot scatterplot matrix for Cluster 2 significant genes
plotClusterSM(2)
# Plot scatterplot matrix for Cluster 3 significant genes
plotClusterSM(3)
# Plot scatterplot matrix for Cluster 4 significant genes
plotClusterSM(4)


