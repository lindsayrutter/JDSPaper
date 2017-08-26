library(pasilla)
library(DESeq2)
library(edgeR)
data("pasillaGenes")
countTable <- counts(pasillaGenes)
countTable <- as.data.frame(countTable)

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


dds <- makeExampleDESeqDataSet(m=6,betaSD=0)
rld <- rlog(dds)
dat0 <- as.data.frame(assay(rld))
p0 <- scatmat(dat0)

dds <- makeExampleDESeqDataSet(m=6,betaSD=1)
rld <- rlog(dds)
dat1 <- as.data.frame(assay(rld))
p1 <- scatmat(dat1)

dds <- makeExampleDESeqDataSet(m=6,betaSD=2)
rld <- rlog(dds)
dat2 <- as.data.frame(assay(rld))
p2 <- scatmat(dat2)

dds <- makeExampleDESeqDataSet(m=6,betaSD=5)
rld <- rlog(dds)
dat5 <- as.data.frame(assay(rld))
p5 <- scatmat(dat5)

dds <- makeExampleDESeqDataSet(m=6,betaSD=5,dispMeanRel = function(x) 0.0000000001/x + .0000000001)
rld <- rlog(dds)
dat5b <- as.data.frame(assay(rld))
p5b <- scatmat(dat5b)

############################################################
dds0 <- makeExampleDESeqDataSet(n=1000,m=6,betaSD=0)
rld0 <- rlog(dds0)
dds1 <- makeExampleDESeqDataSet(n=200,m=6,betaSD=1)
rld1 <- rlog(dds1)
dds2 <- makeExampleDESeqDataSet(n=50,m=6,betaSD=2)
rld2 <- rlog(dds2)
dds3 <- makeExampleDESeqDataSet(n=1000,m=6,betaSD=0.5)
rld3 <- rlog(dds3)
dds4 <- makeExampleDESeqDataSet(n=50,m=6,betaSD=3)
rld4 <- rlog(dds4)

dat0 <- as.data.frame(assay(rld0))
dat1 <- as.data.frame(assay(rld1))
dat2 <- as.data.frame(assay(rld2))
dat3 <- as.data.frame(assay(rld3))
dat4 <- as.data.frame(assay(rld4))
dat <- rbind(dat0,dat1,dat2,dat3,dat4)
p <- scatmat(dat) + xlim(0,13) + ylim(0,13)

############################################################
dds0 <- makeExampleDESeqDataSet(n=1000,m=6,betaSD=0)
rld0 <- rlog(dds0)
dds1 <- makeExampleDESeqDataSet(n=500,m=6,betaSD=1)
rld1 <- rlog(dds1)


dat0 <- as.data.frame(assay(rld0))
dat1 <- as.data.frame(assay(rld1))
dat <- rbind(dat0,dat1)
p <- scatmat(dat)


############################################################
dds0 <- makeExampleDESeqDataSet(n=1000,m=6,betaSD=0)
dds1 <- makeExampleDESeqDataSet(n=500,m=6,betaSD=1)
dds2 <- makeExampleDESeqDataSet(n=500,m=6,betaSD=1.5)
dat <- rbind(assay(dds0), assay(dds1), assay(dds2))
dat <- as.data.frame(dat)
rownames(dat) <- paste0("gene",1:2000)
logDat <- log(dat+1)
p <- scatmat(logDat)
boxplot(logDat)

############################################################

dds1 <- makeExampleDESeqDataSet(n=1000,m=6,betaSD=1)
rld1 <- rlog(dds1)
dat1 <- as.data.frame(assay(rld1))
p <- scatmat(dat1)
