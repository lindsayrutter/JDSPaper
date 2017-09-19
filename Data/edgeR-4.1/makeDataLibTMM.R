# Make dataRaw object
library(readxl)
data <- read_excel("./journal.pone.0009317.s009.xls")
data <- data[3:nrow(data),]
dataSel <- data[,c(4:9)]
colnames(dataSel) <- c("N.8","T.8","N.33","T.33","N.51","T.51")
rownames(dataSel) = as.character(unlist(as.data.frame(data[,1])))
data = as.data.frame(dataSel)
data$N.8 = as.numeric(data$N.8)
data$T.8 = as.numeric(data$T.8)
data$N.33 = as.numeric(data$N.33)
data$T.33 = as.numeric(data$T.33)
data$N.51 = as.numeric(data$N.51)
data$T.51 = as.numeric(data$T.51)
data = data[,c(1,3,5,2,4,6)]
saveRDS(data, "dataRaw.rds")

# Make dataLib object
libSizes <- colSums(data)
dataLib = data.frame('a'=data[,1]/libSizes[1])
for (i in 2:ncol(data)){
  dataLib[,i] = data[,i]/libSizes[i]
}
colnames(dataLib) = colnames(data)
rownames(dataLib) = rownames(data)
setDT(dataLib, keep.rownames = TRUE)[]
colnames(dataLib)[1] = "ID"
dataLib = as.data.frame(dataLib)
saveRDS(dataLib, "dataLib.rds")

# Make dataTMM object
tmmSizes <- calcNormFactors(data, logratioTrim=.3)
dataTMM = data.frame('a'=data[,1]/tmmSizes[1])
for (i in 2:ncol(data)){
  dataTMM[,i] = data[,i]/tmmSizes[i]
}
colnames(dataTMM) = colnames(data)
rownames(dataTMM) = rownames(data)
setDT(dataTMM, keep.rownames = TRUE)[]
colnames(dataTMM)[1] = "ID"
dataTMM = as.data.frame(dataTMM)
saveRDS(dataTMM, "dataTMM.rds")

# Make dataTMM2 object (see if get same TMM normalization factors)
y <- DGEList(counts=data, group= c(rep(1,3),rep(2,3)))
y <- calcNormFactors(y, method="TMM")
identical(tmmSizes, y[[2]]$norm.factors)
identical(data, as.data.frame(y[[1]]))

