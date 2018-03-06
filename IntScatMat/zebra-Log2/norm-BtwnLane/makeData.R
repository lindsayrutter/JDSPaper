library(RUVSeq)
library(zebrafishRNASeq)
library(data.table)

data(zfGenes)

zfGenes <- DGEList(counts=zfGenes)

filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]

filtered.new <- cpm(filtered, TRUE, TRUE)
filtered.norm <- betweenLaneNormalization(filtered.new, which="full", round=FALSE)
filtered = as.data.frame(filtered.norm)

setDT(filtered, keep.rownames = TRUE)[]
colnames(filtered)[1] <- "ID"
filtered$ID <- as.character(filtered$ID)
filtered <- as.data.frame(filtered)
colnames(filtered)[2:7] <- c("C.1", "C.3", "C.5", "T.9", "T.11", "T.13")

zebraData <- filtered

save(zebraData,file="zebraData.Rda")
