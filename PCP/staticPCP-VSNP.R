library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(DESeq2)

outDir = "/Users/lindz/BeeVirusDiet/tblshoot-VSNP/Method1/PCP"

# dat <- read.delim(file="../AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
# colnames(dat) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
# countdata <- dat[ , order(names(dat))]
# countdata <- as.matrix(countdata)
# coldata <- data.frame(row.names = colnames(countdata), virus = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],1,1))), diet = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],2,2))), treatment = unlist(lapply(colnames(countdata), function (x) unlist(strsplit(x, "[.]"))[1])))
# 
# dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
# dds <- DESeq(dds)
# rld <- rlog(dds)

dds <- readRDS("/Users/lindz/BeeVirusDiet/beeDataDDSRLD.rds")[[1]]
rld <- readRDS("/Users/lindz/BeeVirusDiet/beeDataDDSRLD.rds")[[2]]

# Change group1 and group2 as needed
group1 ="NP"
group2 ="VS"

sampleIndex <- which(sapply(colnames(assay(rld)), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1, group2))

bindataSel <- as.data.frame(assay(rld))[, sampleIndex]
setDT(bindataSel, keep.rownames = TRUE)[]
colnames(bindataSel)[1] <- "ID"
bindataSel$ID <- as.factor(bindataSel$ID)
bindataSel <- as.data.frame(bindataSel)

res <- results(dds, contrast=c("treatment",group1,group2))
degIndex <- which(res@listData$padj<0.05) 
pcpDat <- bindataSel[degIndex,]

nVar = ncol(bindataSel)
colNms <- colnames(bindataSel[, c(2:nVar)])
    
boxDat <- bindataSel[, c(1:nVar)] %>% gather(key, val, -c(ID))
colnames(boxDat) <- c("ID", "Sample", "Count")
pcpDat2 <- pcpDat[, c(1:nVar)] %>% gather(key, val, -c(ID))
colnames(pcpDat2) <- c("ID", "Sample", "Count")

p <- ggplot(boxDat, aes(x = Sample, y = Count)) + geom_boxplot() + geom_line(data=pcpDat2, aes(x = Sample, y = Count, group = ID), size = 0.1, alpha = 0.3, color = "red")

jpeg(filename=paste0(outDir, "/", group1, "_", group2, ".jpg"), height=1400, width=1400)
print(p)
dev.off()
