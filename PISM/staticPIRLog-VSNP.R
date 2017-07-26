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
outDir = "/Users/lindz/BeeVirusDiet/tblshoot-VSNP/Method1/PISM"

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
  
  p <- ggplot(data = redPoints, aes_string(x = xChar)) + geom_point(aes_string(y = yChar), size=2, color = "red") + geom_point(data = bluePoints, aes_string(y = yChar), size=1, color = "blue") + geom_point(data = greyPoints, aes_string(y = yChar), size=size, color = "darkgrey") + geom_ribbon(data= pred_interval, aes(ymin = lwr, ymax = upr), fill = "cornflowerblue", alpha = 0.2) + coord_cartesian(xlim = c(minX, maxX), ylim = c(minX, maxX))
  p
}

# Following DESeq2 vignette (using rlog)
dat <- read.delim(file="/Users/lindz/BeeVirusDiet/AllLaneCount.txt",row.names=1,stringsAsFactors = FALSE)
colnames(dat) <- c("NC.1", "NC.2", "NR.1", "VR.1", "NS.1", "VP.1", "NS.2", "VR.2", "NP.1", "VP.2", "VC.1", "NP.2", "VP.3", "NP.3", "VS.1", "VS.2", "VC.2", "NC.3", "VP.4", "NC.4", "NR.2", "VC.3", "VC.4", "NP.4", "VR.3", "NC.5", "VS.3", "NP.5", "VC.5", "VS.4", "NS.3", "VS.5", "VP.5", "NR.3", "NR.4", "VC.6", "NS.4", "NC.6", "NP.6", "VR.4", "NR.5", "NR.6", "NS.5", "VP.6", "NS.6", "VR.5", "VR.6", "VS.6")
countdata <- dat[ , order(names(dat))]
countdata <- as.matrix(countdata)
coldata <- data.frame(row.names = colnames(countdata), virus = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],1,1))), diet = unlist(lapply(colnames(countdata), function (x) substring(unlist(strsplit(x, "[.]"))[1],2,2))), treatment = unlist(lapply(colnames(countdata), function (x) unlist(strsplit(x, "[.]"))[1])))

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)
dds <- DESeq(dds)
rld <- rlog(dds)

# Create prediction interval data frame with upper and lower lines corresponding to sequence covering minimum and maximum of x values in original dataset (don't consider last column because it is not on the x axis in any of the individual scatterplots)
lDat = as.data.frame(assay(rld))
minX <- min(lDat[,c(2:(ncol(lDat)-1))])
maxX <- max(lDat[,c(2:(ncol(lDat)-1))])
newx <- seq(minX - (maxX-minX), maxX + (maxX-minX), by=0.05)

# Change group1 and group2 as needed
group1 ="NP"
group2 ="VS"

sampleIndex <- which(sapply(colnames(assay(rld)), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1, group2))

bindataSel <- as.data.frame(assay(rld))[, sampleIndex]
setDT(bindataSel, keep.rownames = TRUE)[]
colnames(bindataSel)[1] <- "ID"
bindataSel$ID <- as.character(bindataSel$ID)
bindataSel <- as.data.frame(bindataSel)

res <- results(dds, contrast=c("treatment",group1,group2))
degIndex <- which(res@listData$padj<0.05) 
degData <- bindataSel[degIndex,]

maxVal = max(bindataSel[,-1])
minVal = min(bindataSel[,-1])
maxRange = c(minVal, maxVal)

p <- ggpairs(bindataSel[,-1], lower = list(continuous = my_fn))
jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_level99.jpg"), height=1400, width=1400)
print(p)
dev.off()

