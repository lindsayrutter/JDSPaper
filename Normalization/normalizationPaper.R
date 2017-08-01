# 1. R code used in normal distribution test
rawData <- read.table(file.path("input.txt"), header=TRUE)
x <- rawData[,c(1)]
shapiro.test(x)

# 2. R code used in calculating Spearman correlation
inputFile <- "data/uhr.txt"
rawData <- read.table(file.path(inputFile), header=TRUE, sep="\t")
x <- rawData[,c(1)]
a1 <- rawData[,c(2)]
cor.test(x,a1,method="spearman")

#3. R code used in calculating eight non-abundance estimation methods
library(limma)
library(edgeR)
library(DESeq)

inputFile <- " input.txt"
rawData <- read.table(file.path(inputFile), header=TRUE, sep="\t")
geneLength <- rawData[,c(2)]
geneLengthE <- rawData[,c(3)]
geneCount <- rawData[,c(5:11)]

RC <- rowSums(geneCount)

UQ <- apply(geneCount, 1, quantile, 0.75)

Med <- apply(geneCount, 1, median)

f <- calcNormFactors(geneCount, method="TMM")
TMM <- rowSums(geneCount / f)

ef <- estimateSizeFactorsForMatrix(geneCount)
DESeq <- rowSums(geneCount / ef)

MQ <- limma::normalizeQuantiles(geneCount)
Q <- rowSums(MQ)

lengthsInKb <- geneLength/1000
millionMapped <- sum(RC)/1e+06
rpm <- RC/millionMapped
RPKM <- rpm/lengthsInKb

lengthsEInKb <- geneLengthE/1000
ERPKM <- rpm/lengthsEInKb

geneCountNorm <- cbind(RC, UQ, Med, TMM, DESeq, Q, RPKM, ERPKM)

write.table(geneCountNorm, file=" NormResults.txt")

