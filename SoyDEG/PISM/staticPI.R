# Points that are DEG but not outside PI are plotted blue, points that are DEG and outside PI are plotted red, and points that are not DEG but are outside PI are plotted grey.

library(GGally)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(DESeq2)

# Look at VS and NP normalized/filtered to entire dataset
predLevel = 0.999999
size = 0.1
outDir = "/Users/lindz/JDSPaper/PISM"

my_fn <- function(data, mapping, ...){
  xChar = as.character(mapping$x)
  yChar = as.character(mapping$y)
  x = data[,c(xChar)]
  y = data[,c(yChar)]
  m <- lm(y ~ x, data = data)
  mpi <- cbind(data, predict(m, interval = "prediction", level=predLevel))
  # Keep only points that are outside the prediction interval
  plotPoints <- mpi[which(!(mpi[yChar] > mpi$lwr & mpi[yChar] < mpi$upr)),]
  pred_interval <- predict(m, newdata=data.frame(x=newx), interval="prediction", level = predLevel)
  pred_interval <- as.data.frame(pred_interval)
  pred_interval[xChar] = newx
  
  indexBoth = rownames(plotPoints) %in% rownames(degData)
  indexBlue = rownames(degData) %in% rownames(plotPoints)
  redPoints = plotPoints[indexBoth,]
  greyPoints = plotPoints[!indexBoth,] # problem if indexBoth is integer(0)
  bluePoints = degData[!indexBlue,]
  
  p <- ggplot(data = redPoints, aes_string(x = xChar)) + geom_point(aes_string(y = yChar), size=2, color = "red") + geom_point(data = bluePoints, aes_string(y = yChar), size=0.01, alpha=0.5, color = "blue") + geom_point(data = greyPoints, aes_string(y = yChar), size=size, color = "darkgrey") + geom_ribbon(data= pred_interval, aes(ymin = lwr, ymax = upr), fill = "cornflowerblue", alpha = 0.2) + coord_cartesian(xlim = c(minX, maxX), ylim = c(minX, maxX))
  p
}

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
rld <- rlog(dds)

# Create prediction interval data frame with upper and lower lines corresponding to sequence covering minimum and maximum of x values in original dataset (don't consider last column because it is not on the x axis in any of the individual scatterplots)
lDat = as.data.frame(assay(rld))
minX <- min(lDat)
maxX <- max(lDat)
newx <- seq(minX - (maxX-minX), maxX + (maxX-minX), by=0.05)

sampleIndex <- which(sapply(colnames(assay(rld)), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1, group2))

bindataSel <- as.data.frame(assay(rld))[, sampleIndex]
setDT(bindataSel, keep.rownames = TRUE)[]
colnames(bindataSel)[1] <- "ID"
bindataSel$ID <- as.character(bindataSel$ID)
bindataSel <- as.data.frame(bindataSel)

res <- results(dds, contrast=c("treatment",group1,group2))
degIndex <- which(res@listData$padj<0.05) 
degData <- bindataSel[degIndex,]

maxRange = c(minX, maxX)

p <- ggpairs(bindataSel[,-1], lower = list(continuous = my_fn))
jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_level99.jpg"), height=1400, width=1400)
print(p)
dev.off()

