source("functions.R")

# ------------------------------------
# table of counts from Marioni et al.
# ------------------------------------
load("LK_data.RData")
D <- as.matrix(MA.subsetA$M)
g <- as.character(MA.subsetA$genes$EnsemblGeneID)
o <- order(gsub("R[1-2]L[1-8]","",colnames(D)))

# ------------------------------------
# download symbols and RefSeq ids from (ERROR)
# ------------------------------------
if(!file.exists("bm.Rdata")) {
  library(biomaRt)
  mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
  
  # get refseq identifiers from biomaRt
  bm <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','refseq_dna'), filters="ensembl_gene_id", value=g, mart=mart)
  save(bm, file="bm.Rdata")
} else {
  load("bm.Rdata")
}

# ------------------------------------
# read 575 housekeeping genes from file
# ------------------------------------
hk <- read.table("human_housekeeping.txt", header=FALSE, sep=" ", as.is=TRUE, stringsAsFactors=FALSE)$V1
hk1 <- gsub("*","",hk, fixed=TRUE)

# ------------------------------------
# match housekeeping genes to our table of data (ERROR)
# ------------------------------------
m <- match(hk1, bm$refseq_dna)
m <- m[!is.na(m)]
ensHk <- unique(bm$ensembl_gene_id[m])
m <- match(ensHk, g)


# ------------------------------------
# simulate data from empirical distribution of counts
# ------------------------------------

libSizes <- colSums(D)

foldDiff <- 2
pUp <- .8
pDifferential=0.1

xx <- generateDataset(commonTags=2e4, uniqueTags=c(3000,100), foldDifference=foldDiff, pUp=pUp, pDifferential=pDifferential, empiricalDist=D[,1], libLimits=rep(1e6,2))

ci <- xx$commonInd
di <- xx$differentialInd
nDiff <- length(di)


library(edgeR)
fW <- calcNormFactors(xx$DATA)[2]

library(statmod)
# ------------------------------------
# calculate Fisher test using library size adjustment
# ------------------------------------
s <- sage.test(xx$DATA[,1], xx$DATA[,2], n1=sum(xx$DATA[,1]), n2=sum(xx$DATA[,2]))
o <- order(s)
w <- o %in% di

# ------------------------------------
# calculate Fisher test after adjusting
# ------------------------------------
sM <- sage.test(xx$DATA[,1], xx$DATA[,2], n1=sum(xx$DATA[,1])/sqrt(fW), n2=sum(xx$DATA[,2])*sqrt(fW))
oM <- order(sM)
wM <- oM %in% di

################# Lindsay: Applied this for TMM to real data ################# 
f <- calcNormFactors(D, logratioTrim=.3) 

ff1 <- f[3]
ff2r <- f[2]

cols <- rep("darkgray",nrow(xx$DATA))
cols[intersect(di,ci)] <- "blue"
cols[setdiff(di,ci)] <- "orange"


ma1 <- maPlot(D[,1], D[,2], normalize=TRUE, pch=19, cex=.5,ylim=c(-8,8), allCol="darkgray", xlab=expression( A == log[2] (sqrt(Liver %.% Kidney))  ), ylab=expression(M == log[2](Liver/Kidney)))
dev.off()

layoutMatrix <- matrix( c(1,2,4,4,3,3,5,5), nc=2 )


# ------------------------------------
# Figure 1
# ------------------------------------
pdf("GB_fig1.pdf",8,4)
nf <- layout( layoutMatrix[1:2,] )
par(mai=c(0.6732,0.639,0.1,0.05))

################# Lindsay: Here is where do normalization by library size ################# 
hist( log2( (D[,1]/libSizes[1]) / (D[,3]/libSizes[3]) ), 20, main="", xlab=expression(log[2](Kidney1/N[K1])-log[2](Kidney2/N[K2])), prob=TRUE ,xlim=c(-7,7))
abline(v=log2(ff1), lwd=3, col="red")
mtext("A", side=3, adj=-.16, padj=1, cex=1.5)

hist( -log2( (D[,1]/libSizes[1]) / (D[,2]/libSizes[2]) ), 50, main="", xlab=expression(log[2](Liver/N[L])-log[2](Kidney/N[K])), prob=TRUE,xlim=c(-7,7), ylim=c(0,.45) )
lines(density(ma1$M[m]), col="green",lwd=2)
abline(v=0, lty=3, col="black", lwd=3)
mtext("B", side=3, adj=-.16, padj=1, cex=1.5)
abline(v=log2(ff2r), lwd=3, col="red")

ma1 <- maPlot(D[,1], D[,2], normalize=TRUE, pch=19, cex=.5,ylim=c(-8,8),allCol="darkgray", xlab=expression( A == log[2] (sqrt(Liver/N[L] %.% Kidney/N[K]))  ), ylab=expression(M == log[2](Liver/N[L])-log[2](Kidney/N[K])))
grid(col="black")
points(ma1$A[m], ma1$M[m], pch=19, col="green", cex=.4)
abline(h=median(ma1$M[m]), col="green",lwd=1.5)
abline(h=log2(ff2r), col="red",lwd=1.5)
legend("bottomright",legend=c("housekeeping genes","unique to a sample"),col=c("green","orange"),pch=19,cex=0.9,bg="white")
arrows( -9, 8, -10.5, 7, length=.1, lwd=4 )
mtext("C", side=3, adj=-.16, padj=1, cex=1.5)
dev.off()


# ------------------------------------
# Figure 2
# ------------------------------------
pdf("GB_fig2.pdf",h=4,7)
par(mfrow=c(1,2))
par(mai=c(0.782,0.789,0.1,0.1))

maPlot(xx$DATA[,1], xx$DATA[,2], normalize=TRUE, pch=".", cex=2.5, col=cols, allCol="darkgray", ylim=c(-8,8))
grid(col="black")
abline(h=log2(fW), col="red")
legend("bottomright",legend=c("2 fold DE","unique to a sample"),col=c("blue","orange"),pch=".",pt.cex=3,cex=0.9,bg="white")
mtext("A", side=3, adj=-.25, padj=1, cex=1.5)

plot(1:nDiff, cumsum(!w[1:nDiff]), xlab="Number of Genes Selected", ylab="Number of False Discoveries", lwd=4, type="l" )
lines(1:nDiff, cumsum(!wM[1:nDiff]), lwd=4, type="l", col="red" )
grid(col="black")
legend("topleft",c("Fisher Test (total reads)","Fisher Test (TMM)"),col=c("black","red"),lwd=4, bg="white",cex=0.8)
mtext("B", side=3, adj=-.25, padj=1, cex=1.5)

dev.off()