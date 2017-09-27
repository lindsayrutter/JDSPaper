library(RUVSeq)
library(zebrafishRNASeq)
library(data.table)

data(zfGenes)

zfGenes <- DGEList(counts=zfGenes)

filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]

filtered <- calcNormFactors(filtered, method = "TMM")
rld <- rlog(filtered[[1]])
rld <- as.data.frame(rld)

setDT(rld, keep.rownames = TRUE)[]
colnames(rld)[1] <- "ID"
rld$ID <- as.character(rld$ID)
rld <- as.data.frame(rld)
colnames(rld)[2:7] <- c("C.1", "C.3", "C.5", "T.9", "T.11", "T.13")

zebraData <- rld

save(zebraData,file="zebraData.Rda")
