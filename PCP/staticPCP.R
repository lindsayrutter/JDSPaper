library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(DESeq2)

outDir = "/Users/lindz/JDSPaper/PCP"
load("../All_leaves040615.rda")
ct <- countTable
ct2 <- assays(ct)[[1]]
rownames(ct2) <- ct@rowRanges@elementMetadata@listData$ID
colnames(ct2) <- unlist(strsplit(colnames(ct2), "\\."))[seq(1, 17*3, 3)]
leaves.all <- ct2
bindataSel = as.data.frame(leaves.all[,c("ML08R","ML14R","ML22R","ML11R","ML27R","ML33R")])
colnames(bindataSel) <- c("N.1","N.2","N.3","P.1","P.2","P.3")
bindataSel <- as.data.frame(bindataSel)

coldata <- data.frame(row.names = colnames(bindataSel), treatment = unlist(lapply(colnames(bindataSel), function (x) unlist(strsplit(x, "[.]"))[1])))
dds <- DESeqDataSetFromMatrix(countData = bindataSel, colData = coldata, design = ~ treatment)
dds <- DESeq(dds)
group1 <- levels(coldata$treatment)[1]
group2 <- levels(coldata$treatment)[2]

res <- results(dds, contrast=c("treatment",group1,group2))
degIndex <- which(res@listData$padj<0.0001)
# Or take top 5 most differentially expressed
degIndex <- which(res@listData$padj %in% sort(res@listData$padj)[1:15])
logBin <- as.data.frame(assay(rlog(dds)))

setDT(logBin, keep.rownames = TRUE)[]
colnames(logBin)[1] <- "ID"
logBin$ID <- as.character(logBin$ID)
logBin <- as.data.frame(logBin)
pcpDat <- logBin[degIndex,]
    
boxDat <- logBin %>% gather(key, val, -c(ID))
colnames(boxDat) <- c("ID", "Sample", "Count")
pcpDat2 <- pcpDat %>% gather(key, val, -c(ID))
colnames(pcpDat2) <- c("ID", "Sample", "Count")

p <- ggplot(boxDat, aes(x = Sample, y = Count)) + geom_boxplot() + geom_line(data=pcpDat2, aes(x = Sample, y = Count, group = ID), size = 0.1, color = "red") #alpha = 0.3

jpeg(filename=paste0(outDir, "/", group1, "_", group2, ".jpg"), height=1400, width=1400)
print(p)
dev.off()
