library(pasilla)
library(DESeq2)
library(edgeR)
data("pasillaGenes")
countTable <- counts(pasillaGenes)
countTable <- as.data.frame(countTable)

colnames(countTable) <- c("T.1","T.2","T.3","U.1","U.2","U.3","U.4")
setDT(countTable, keep.rownames = TRUE)[]
colnames(countTable)[1] <- "ID"
countTable <- as.data.frame(countTable)
countTable[,-1] <- log(countTable[,-1]+1)
plotScatterStatic(countTable)

y <- DGEList(counts=countTable)
allGroups <- c(rep("treated",3), rep("untreated",3))
plotMDS(y, col = c("blue","blue","blue","blue","red","red","red","red")[factor(allGroups)], cex=0.6)

plotMDS(y, col = c("blue","blue","blue","blue","red","red","red","red")[factor(allGroups)], cex=0.6,
xlim=c(-1.3, 1.3), ylim=c(-1.3, 1.3))

plotMDS(y, col = c("blue","blue","blue","blue","red","red","red","red")[factor(allGroups)], cex=0.6, xlim=c(-1.3, 1.3), ylim=c(-1.3, 1.3)) + theme(aspect.ratio=1)

yd <- as.data.frame(y[[1]])
ly <- log(yd+1)
scatmat(ly)

rldY <- rlog(y)

