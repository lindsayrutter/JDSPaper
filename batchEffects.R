# From https://support.bioconductor.org/p/93333/
library(tidyverse)
library(genefilter)
# This cannot be used for scatterplot Matrices

n <- 10000
m <- 20
x <- matrix(rnorm(n*m), nrow = n, ncol = m)
fac <- factor(c(rep(0, 10), rep(1, 10)))
# Ttests for rows of matrix
rt1 <- rowttests(x, fac)
# Data looks okay
qplot(rt1$p.value, fill = I("tan3")) #batchEffectImage1.png

# Add artificial batch effect by adding 5 to groups 6-10 and 11-15
x[, 6:15] <- x[, 6:15]+0.5
rt2 <- rowttests(x, fac)
qplot(rt2$p.value, fill = I("coral3")) #batchEffectImage2.png

myData = as.data.frame(x)
myData = myData[,c(4:9)]
colnames(myData) = c("N.1","N.2","N.3","B.1","B.2","B.3")
myData = abs(myData)
setDT(myData, keep.rownames = TRUE)[]
colnames(myData)[1] <- "ID"
myData = as.data.frame(myData)
plotScatterStatic(myData)

myData[,-1] = log(myData[,-1]+1)
plotScatterStatic(myData)

