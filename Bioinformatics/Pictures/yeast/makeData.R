library(devtools)
install_github("drisso/yeastRNASeqRisso2011")
library(yeastRNASeqRisso2011)
#library(yeastRNASeq)
library(EDASeq)

####################################################################
# I begin with trying to copy the same methodology seen in the EDASeq vignette. However, this causes an error at withinLaneNormalization

load("~/JDSPaper/Data/risso-Yeast/data/geneLevelCounts.rda")
load("~/JDSPaper/Data/risso-Yeast/data/laneInfo.rda")
load("~/JDSPaper/Data/risso-Yeast/data/geneInfo.rda")
means <- rowMeans(geneLevelCounts)
filter <- means >= 10
table(filter)
geneLevelCounts <- geneLevelCounts[filter,]
data <- newSeqExpressionSet(counts = as.matrix(geneLevelCounts), featureData = geneInfo[rownames(geneLevelCounts), ], phenoData = laneInfo)

# Plot checks
colors <- as.numeric(pData(data)[, 2]) + 1
boxplot(data, col=colors)
meanVarPlot(data[, 1:8], log=TRUE)
biasPlot(data[,1:8], "gc", log=TRUE, ylim=c(0, 8), col=1) # Error

# There are 4 within-normalization types, and 3 between-normalization types. Here, we do one of each
dataWithin <- withinLaneNormalization(data, "gc", which="full", offset=FALSE)

####################################################################
# This vignette methodology results in an error. So, I start over using some functions and data from the help function of withinLaneNormalization along with the vignette
rm(list=ls())

load("~/JDSPaper/Data/risso-Yeast/data/geneLevelCounts.rda")
load("~/JDSPaper/Data/risso-Yeast/data/laneInfo.rda")
load("~/JDSPaper/Data/risso-Yeast/data/geneInfo.rda")
data(yeastGC)

#Need to change column name otherwise get an error at future command (counts <- as(dataNorm,"CountDataSet"))
colnames(laneInfo)[2] <- "conditions"

means <- rowMeans(geneLevelCounts)
filter <- means >= 10
table(filter)
geneLevelCounts <- geneLevelCounts[filter,]

sub <- intersect(rownames(geneLevelCounts), names(yeastGC))
mat <- as.matrix(geneLevelCounts[sub, ])
data <- newSeqExpressionSet(mat, phenoData=laneInfo, featureData=AnnotatedDataFrame(data.frame(gc=yeastGC[sub])))

colors <- as.numeric(pData(data)[, 2]) + 1
boxplot(data, col=colors)
meanVarPlot(data[, 1:8], log=TRUE)
biasPlot(data[,1:8], "gc", log=TRUE, ylim=c(0, 8), col=1) # Error

dataWithin <- withinLaneNormalization(data, "gc", which="full", offset=FALSE)
dataNorm <- betweenLaneNormalization(dataWithin, which="median")

biasPlot(dataNorm[,1:8], "gc", log=TRUE, ylim=c(0, 8), col=1)
boxplot(dataNorm, col=colors)
meanVarPlot(dataNorm[, 1:8], log=TRUE)

####################################################################
# FYI: How to transform the EDASeq normalized data to be used in edgeR/DESeq for differential analysis
library(DESeq)
counts <- as(dataNorm,"CountDataSet")
str(as.data.frame(counts@assayData$counts))
y <- DGEList(counts=counts@assayData$counts)
