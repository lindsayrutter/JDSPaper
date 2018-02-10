library(rtracklayer)
library(Rsamtools)
library(grid)
library(GenomicAlignments)
library(ggplot2)
library(GGally)
library(edgeR)
library(stringr)
library(EDASeq)
library(dplyr)
library(matrixStats)
library(gridExtra)
library(reshape2)
library(scales)
library(bigPint)

data("soybean_ir")
data("soybean_ir_metrics")
data <- soybean_ir
metrics <- soybean_ir_metrics[["N_P"]]

RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

# Make sure each gene has at least one count in at least half of the six samples
filterLow = which(rowSums(data[,-1])<=ncol(data[,-1])/2)
filt1 <- data[filterLow,]
rownames_filt1 <- filt1$ID
filt1 <- filt1[,-1]
filt1 = mutate(filt1, mean = (N.1+N.2+N.3+P.1+P.2+P.3)/6, stdev = RowSD(cbind(N.1,N.2,N.3,P.1,P.2,P.3)))
rownames(filt1) <- rownames_filt1

data <- data[-filterLow,]
data_Rownames <- data$ID
data = data[,-1]
rownames(data) <- data_Rownames
#Normalize and log
cpm.data.new <- cpm(data, TRUE, TRUE)
# Normalize for sequencing depth and other distributional differences between lanes
data <- betweenLaneNormalization(cpm.data.new, which="full", round=FALSE)
data = as.data.frame(data)
# Add mean and standard deviation for each row/gene
data = mutate(data, mean = (N.1+N.2+N.3+P.1+P.2+P.3)/6, stdev = RowSD(cbind(N.1,N.2,N.3,P.1,P.2,P.3)))
rownames(data)=data_Rownames
data$ID <- data_Rownames
# Remove the genes with lowest quartile of mean and standard deviation
qT = as.numeric(summary(data$mean)["1st Qu."])
dataq = subset(data,mean>qT)
qTs = as.numeric(summary(dataq$stdev)["1st Qu."])
dataq = subset(dataq,stdev>qTs)
filt = subset(data,mean<=qT|stdev<=qTs)
filt <- rbind(filt[,-9], filt1)
filt$ID <- rownames(filt)

# Apply Loess model and further filter low gene counts
model = loess(mean ~ stdev, data=dataq)
dataqp = dataq[which(sign(model$residuals) == 1),]
dataqn = dataq[which(sign(model$residuals) == -1),]
dataqp = dataqp[,1:6]

#Scale filter data
filt = filt[,1:6]
filt = rbind(filt,dataqn[,1:6])

dataqps <- t(apply(as.matrix(dataqp[,1:6]), 1, scale))
filts <- t(apply(as.matrix(filt[,1:6]), 1, scale))
dataqps <- as.data.frame(dataqps)
colnames(dataqps) <- colnames(dataqp[,1:6])
dataqps$ID <- rownames(dataqps)
filts <- as.data.frame(filts)
colnames(filts) <- colnames(filt[,1:6])
filts$ID <- rownames(filts)
# Indices of the 9760 NAN rows. They had stdev=0 in the filt data
nID <- which(is.nan(filts$N.1))
# Set these filtered values that have all same values for samples to 0
filts[nID,1:6] <- 0

# Comine the filtered and remaining data
fulls <- rbind(dataqps, filts)
filteredID <- filts$ID
