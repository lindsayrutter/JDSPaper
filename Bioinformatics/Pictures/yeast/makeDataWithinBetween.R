library(devtools)
install_github("drisso/yeastRNASeqRisso2011")
library(yeastRNASeqRisso2011)
library(EDASeq)
library(DESeq)

load("~/JDSPaper/Data/risso-Yeast/data/geneLevelCounts.rda")
load("~/JDSPaper/Data/risso-Yeast/data/laneInfo.rda")
load("~/JDSPaper/Data/risso-Yeast/data/geneInfo.rda")
data(yeastGC)

colnames(laneInfo)[2] <- "conditions"

means <- rowMeans(geneLevelCounts)
filter <- means >= 10
table(filter)
geneLevelCounts <- geneLevelCounts[filter,]

sub <- intersect(rownames(geneLevelCounts), names(yeastGC))
mat <- as.matrix(geneLevelCounts[sub, ])
data <- newSeqExpressionSet(mat, phenoData=laneInfo, featureData=AnnotatedDataFrame(data.frame(gc=yeastGC[sub])))

dataWithin <- withinLaneNormalization(data, "gc", which="full", offset=FALSE)
dataLoess <- withinLaneNormalization(data, "gc", which="loess", offset=FALSE)
dataNorm <- betweenLaneNormalization(dataWithin, which="median")
dataBetweenOnly <- betweenLaneNormalization(data, which="median")
dataLoessBetween <- betweenLaneNormalization(dataLoess, which="median")

counts <- as(dataBetweenOnly,"CountDataSet")
dataBetweenOnly <- as.data.frame(counts@assayData$counts)
counts <- as(dataNorm,"CountDataSet")
dataBetween <- as.data.frame(counts@assayData$counts)
counts2 <- as(dataWithin,"CountDataSet")
dataWithin <- as.data.frame(counts2@assayData$counts)
dataRaw <- as.data.frame(mat)
counts <- as(dataLoessBetween,"CountDataSet")
dataLoessBetween <- as.data.frame(counts@assayData$counts)

saveRDS(dataRaw, "dataRaw.rds")
saveRDS(dataBetween, "dataBetween.rds")
saveRDS(dataWithin, "dataWithin.rds")
saveRDS(dataBetweenOnly, "dataBetweenOnly.rds")
saveRDS(dataLoessBetween, "dataLoessBetween.rds")
