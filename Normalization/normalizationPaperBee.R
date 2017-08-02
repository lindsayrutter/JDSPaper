library(edgeR)
library(DESeq2)

load("bee-Log2/beeData.Rda")
dat <- beeData
rm(beeData)

RC <- rowSums(dat[,2:49])
UQ <- apply(dat[,2:49], 1, quantile, 0.75)
Med <- apply(dat[,2:49], 1, median)

# edgeR::calcNormFactors
f <- calcNormFactors(dat[,2:49], method="TMM")
TMM <- rowSums(dat[,2:49] / f)

# DESeq2::estimateSizeFactorsForMatrix
ef <- estimateSizeFactorsForMatrix(dat[,2:49])

DESeq <- rowSums(dat[,2:49] / ef)

MQ <- limma::normalizeQuantiles(dat[,2:49])
Q <- rowSums(MQ)

lengthsInKb <- geneLength/1000
millionMapped <- sum(RC)/1e+06
rpm <- RC/millionMapped
RPKM <- rpm/lengthsInKb
