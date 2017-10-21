library(limma)
source("functions.R")
library(readxl)
g <- read_excel("~/JDSPaper/Data/robinson-Cloonan/Grimmond_lengths.xls")

# ------------------------------------
# generate data with replicates from joint 
# distribution of lengths/counts
# ------------------------------------

pUp <- .8
foldDiff <- 4

xx <- generateDataset2(commonTags=2e4, uniqueTags=c(2200,0), group=c(1,2), foldDifference=2, pUp=pUp, pDifferential=0.05, empiricalDist=g$EB, libLimits=c(.6e6,1e6), lengthDist=g$transcript_length)
k <- which(rowSums(xx$DATA) > 0)
xx <- takeSubset(xx, k)

# ------------------------------------
# run Marioni analysis
# ------------------------------------
gg <- xx$group
gl <- levels(xx$group)
pm <- Poisson.model.new(xx$DATA, which(gg==gl[1]), which(gg==gl[2]), calcFactor=FALSE)
pmAdj <- Poisson.model.new(xx$DATA, which(gg==gl[1]), which(gg==gl[2]), calcFactor=TRUE)


# ------------------------------------
# run edgeR-like analysis
# ------------------------------------
fW <- calcNormFactors(xx$DATA)

lambda <- rowSums(xx$DATA)/sum(xx$DATA)
expMean <- outer(lambda, colSums(xx$DATA))
cs <- colSums(xx$DATA)
effM <- cs*fW
expMeanAdj <- outer(lambda, effM)

exactP <- exactTestPoisson(xx$DATA, group1Ind=1:2, group2Ind=3:4, expMean, verbose = TRUE)
exactPadj <- exactTestPoisson(xx$DATA, group1Ind=1:2, group2Ind=3:4, expMeanAdj, verbose = TRUE)

# ------------------------------------
# run Cloonan analysis
# ------------------------------------
library(limma)
d <- log2(xx$DATA+1)/outer(xx$length,rep(1,ncol(xx$DATA)))
d <- normalizeQuantiles(d)

mm <- model.matrix(~xx$group)
f <- lmFit(d, mm)
f <- eBayes(f)
cloonanpval <- f$p.value[,2]

ci <- xx$commonInd
di <- intersect(ci,xx$differentialInd)

save.image("simulations2.Rdata")

# false discoveries amongst the common genes
pdf("GB_fig3.pdf", w=6, h=6)

fdPlot( pm$stats$pval[ci], di, lwd=1 )
fdPlot( pmAdj$stats$pval[ci], di, add=TRUE, lty=2, lwd=1 )

fdPlot( exactP[ci], di, add=TRUE, col="blue", lwd=1 )
fdPlot( exactPadj[ci], di, add=TRUE, lty=2, col="blue", lwd=1 )

fdPlot( cloonanpval[ci], di, add=TRUE, col="red", lwd=1 )

legend("topleft", c("Cloonan","Poisson-LR","Poisson-LR-TMMnormalized","Poisson-Exact",
                    "Poisson-Exact-TMMnormalized"),
       col=c("red","black","black","blue","blue","orange","orange"),
       lty=c(1,1,2,1,2,1,2), lwd=1, cex=.7)

dev.off()


# ------------------------------------
# supplementary figure
# ------------------------------------
cols <- rep("black",nrow(xx$DATA))
cols[setdiff(xx$differentialInd,xx$commonInd)] <- "orange"
cols[intersect(xx$differentialInd,xx$commonInd)] <- "blue"

types <- c("common","common+differential","unique")
names(types) <- c("black","blue","orange")

fW <- calcNormFactors(xx$DATA)

png("figS3.maplots.png", w=1000, h=500)
par(mfrow=c(1,3))
maPlot( xx$DATA[,1], xx$DATA[,2], normalize=TRUE, ylim=c(-5,5), col=cols, pch=19, cex=.3 ); grid(col="blue"); abline(h=log2(fW[2]), lwd=4, col="red")
maPlot( xx$DATA[,1], xx$DATA[,3], normalize=TRUE, ylim=c(-5,5), col=cols, pch=19, cex=.3 ); grid(col="blue"); abline(h=log2(fW[3]), lwd=4, col="red")
maPlot( xx$DATA[,1], xx$DATA[,4], normalize=TRUE, ylim=c(-5,5), col=cols, pch=19, cex=.3 ); grid(col="blue"); abline(h=log2(fW[4]), lwd=4, col="red")
dev.off()