library(data.table)
load("LK_data.RData")
data = as.data.frame(MA.subsetA$M)
rownames(data) = as.character(MA.subsetA$genes$EnsemblGeneID)
data = as.data.frame(data)
data = data[,c(1,3,6,8,10,2,4,5,7,9)]

libSizes <- colSums(data)
dataLib = data.frame('a'=data[,1]/libSizes[1])
for (i in 2:ncol(data)){
  dataLib[,i] = data[,i]/libSizes[i]
}
colnames(dataLib) = colnames(data)
rownames(dataLib) = rownames(data)
setDT(dataLib, keep.rownames = TRUE)[]
colnames(dataLib) = c("ID","K.R1L1","K.R1L3","K.R1L7","K.R2L2","K.R2L6","L.R1L2","L.R1L4","L.R1L6","L.R1L8","L.R2L3")
dataLib = as.data.frame(dataLib)

tmmSizes <- calcNormFactors(data, logratioTrim=.3)
dataTMM = data.frame('a'=data[,1]/tmmSizes[1])
for (i in 2:ncol(data)){
  dataTMM[,i] = data[,i]/tmmSizes[i]
}
colnames(dataTMM) = colnames(data)
rownames(dataTMM) = rownames(data)
setDT(dataTMM, keep.rownames = TRUE)[]
colnames(dataTMM) = c("ID","K.R1L1","K.R1L3","K.R1L7","K.R2L2","K.R2L6","L.R1L2","L.R1L4","L.R1L6","L.R1L8","L.R2L3")
dataTMM = as.data.frame(dataTMM)

saveRDS(dataLib, "dataLib.rds")
saveRDS(dataTMM, "dataTMM.rds")
