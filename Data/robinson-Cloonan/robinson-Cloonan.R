# This is the data from reference 12 of TMM Robinson paper

library(data.table)


load("LK_data.RData")
data = as.data.frame(MA.subsetA$M)
rownames(data) = as.character(MA.subsetA$genes$EnsemblGeneID)
setDT(data, keep.rownames = TRUE)[]
colnames(data) = c("ID","K.R1L1","L.R1L2","K.R1L3","L.R1L4","L.R1L6","K.R1L7","L.R1L8","K.R2L2","L.R2L3","K.R2L6")
data = as.data.frame(data)
data = data[,c(1,2,4,7,9,11,3,5,6,8,10)]

# Obtain 2K 2L
data <- data[,c(1:4,7:9)]
data[,c(2:ncol(data))] = log(data[,c(2:ncol(data))]+1)
plotScatterStatic(data)



