library(edgeR)
library(EDASeq)

# https://rdrr.io/cran/ssizeRNA/man/sim.counts.html
data <- sim.counts(nGenes = 10000, pi0 = 0.8, m=3, mu=125, disp=1, fc=2, up = 0.5, replace = TRUE)
data <- as.data.frame(data$counts)
colnames(data) <- c("A.1","A.2","A.3","B.1","B.2", "B.3")

# Normalize and log
cpm.data.new <- cpm(data, TRUE, TRUE)
# Normalize for sequencing depth and other distributional differences between lanes
cpm.data.norm <- betweenLaneNormalization(cpm.data.new, which="full", round=FALSE)
data = as.data.frame(cpm.data.norm)

# Function to calculate row standard deviation
RowSD = function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}

# Add mean and standard deviation for each row/gene
data = mutate(data, mean = (A.1+A.2+A.3+B.1+B.2+B.3)/6, stdev = RowSD(cbind(A.1,A.2,A.3,B.1,B.2,B.3)))
rownames(data) <- paste0("Gene", 1:nrow(data))

# Remove the genes with lowest quartile of mean and standard deviation
q1T = as.numeric(summary(data$mean)["1st Qu."])
dataq1 = subset(data,mean>q1T)
q1Ts = as.numeric(summary(dataq1$stdev)["1st Qu."])
dataq1 = subset(dataq1,stdev>q1Ts)
# Keep the filtered genes
filt = subset(data,mean<=q1T|stdev<=q1Ts)

# Remove genes with Loess residuals of -1
model = loess(mean ~ stdev, data=dataq1)
dataq1p = dataq1[which(sign(model$residuals) == 1),]
dataq1n = dataq1[which(sign(model$residuals) == -1),]
dataq1p = dataq1p[,1:6]
# Scale each row to have mean=0 and stdev=1
dataq1s = t(apply(as.matrix(dataq1p), 1, scale))
colnames(dataq1s)=colnames(dataq1p)

# Scale each row of filtered genes to have mean=0 and stdev=1
filt = filt[,1:6]
dataq1n = dataq1n[,1:6]
filt = rbind(filt,dataq1n)
filts = t(apply(as.matrix(filt), 1, scale))

dataq1s$ID = rownames(dataq1s)
boxDat <- melt(dataq1s, id.vars="ID")
colnames(boxDat) <- c("ID","Sample","Value")
ggplot(boxDat, aes(x=Sample, y=Value)) + geom_boxplot()




# Fit dendogram of scale genes
dendo = dataq1s
rownames(dendo) = NULL
d = dist(as.matrix(dendo))
hc = hclust(d, method="ward.D")