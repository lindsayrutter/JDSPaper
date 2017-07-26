library(plotly)
library(GGally)
library(hexbin)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(DESeq2)

# Look at VS and NP normalized/filtered to entire dataset

outDir = "/Users/lindz/JDSPaper/HexSM"

my_fn <- function(data, mapping, ...){
  xChar = as.character(mapping$x)
  yChar = as.character(mapping$y)
  x = data[,c(xChar)]
  y = data[,c(yChar)]
  h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
  hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
  attr(hexdf, "cID") <- h@cID
  p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer)) + geom_point(data = degData, aes_string(x=xChar, y=yChar), inherit.aes = FALSE, color = "orange", size = 0.01, shape = ".")
  p
}

load("../All_leaves040615.rda")
ct <- countTable
ct2 <- assays(ct)[[1]]
rownames(ct2) <- ct@rowRanges@elementMetadata@listData$ID
colnames(ct2) <- unlist(strsplit(colnames(ct2), "\\."))[seq(1, 17*3, 3)]
leaves.all <- ct2
bindataSel = as.data.frame(leaves.all[,c("ML08R","ML14R","ML22R","ML11R","ML27R","ML33R")])

# setDT(bindataSel, keep.rownames = TRUE)[]
# colnames(bindataSel)[1] <- "ID"
# bindataSel$ID <- as.character(bindataSel$ID)
# bindataSel <- as.data.frame(bindataSel)
colnames(bindataSel) <- c("N.1","N.2","N.3","P.1","P.2","P.3")
bindataSel <- as.data.frame(bindataSel)

coldata <- data.frame(row.names = colnames(bindataSel), treatment = unlist(lapply(colnames(bindataSel), function (x) unlist(strsplit(x, "[.]"))[1])))
dds <- DESeqDataSetFromMatrix(countData = bindataSel, colData = coldata, design = ~ treatment)
dds <- DESeq(dds)
group1 <- levels(coldata$treatment)[1]
group2 <- levels(coldata$treatment)[2]

res <- results(dds, contrast=c("treatment",group1,group2))
degIndex <- which(res@listData$padj<0.0001) 
#degData <- bindataSel[degIndex,]

logBin <- as.data.frame(assay(rlog(dds)))
degData <- logBin[degIndex,]

maxVal = max(logBin[,-1])
minVal = min(logBin[,-1])
maxRange = c(minVal, maxVal)
xbins=10
buffer = maxRange[2]/xbins
p <- ggpairs(as.data.frame(assay(rlog(dds))), lower = list(continuous = my_fn))
jpeg(filename=paste0(outDir, "/", group1, "_", group2, ".jpg"), height=1400, width=1400)
print(p)
dev.off()

