library(RUVSeq)
library(zebrafishRNASeq)
library(data.table)

data(zfGenes)

zfGenes <- DGEList(counts=zfGenes)

filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]

filtered.norm <- calcNormFactors(filtered, method = "none")
filtered.normlog <- cpm(filtered.norm, log=TRUE)
filtered.normlog = as.data.frame(filtered.normlog)

setDT(filtered.normlog, keep.rownames = TRUE)[]
colnames(filtered.normlog)[1] <- "ID"
filtered.normlog$ID <- as.character(filtered.normlog$ID)
filtered.normlog <- as.data.frame(filtered.normlog)
colnames(filtered.normlog)[2:7] <- c("C.1", "C.3", "C.5", "T.9", "T.11", "T.13")

zebraData <- filtered.normlog

save(zebraData,file="zebraData.Rda")
