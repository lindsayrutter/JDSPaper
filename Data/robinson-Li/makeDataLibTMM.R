library(data.table)
library(readxl)
data1 <- read_excel("~/JDSPaper/Data/robinson-Cloonan/Grimmond_lengths.xls")

s1 <- read_excel("~/JDSPaper/Data/robinson-Li/li_mapped/s1.xls")
s2 <- read_excel("~/JDSPaper/Data/robinson-Li/li_mapped/s2.xls")
s3 <- read_excel("~/JDSPaper/Data/robinson-Li/li_mapped/s3.xls")
s4 <- read_excel("~/JDSPaper/Data/robinson-Li/li_mapped/s4.xls")
s5 <- read_excel("~/JDSPaper/Data/robinson-Li/li_mapped/s5.xls")
s6 <- read_excel("~/JDSPaper/Data/robinson-Li/li_mapped/s6.xls")
s8 <- read_excel("~/JDSPaper/Data/robinson-Li/li_mapped/s8.xls")

colnames(s1) <- c("ID", "S.1")
colnames(s2) <- c("ID", "S.2")
colnames(s3) <- c("ID", "S.3")
colnames(s4) <- c("ID", "S.4")
colnames(s5) <- c("ID", "S.5")
colnames(s6) <- c("ID", "S.6")
colnames(s8) <- c("ID", "S.8")

data <- as.data.frame(cbind(s1$S.1, s2$S.2, s3$S.3, s4$S.4, s5$S.5, s6$S.6, s8$S.8))
data = as.data.frame(data)
colnames(data) = c("S.1","S.2","S.3","S.4","S.5","S.6","S.8")

# Make dataLib object
libSizes <- colSums(data)
dataLib = data.frame('a'=data[,1]/libSizes[1])
for (i in 2:ncol(data)){
  dataLib[,i] = data[,i]/libSizes[i]
}
colnames(dataLib) = colnames(data)
rownames(dataLib) = rownames(data)
setDT(dataLib, keep.rownames = TRUE)[]
colnames(dataLib) = c("ID","S.1","S.2","S.3","S.4","S.5","S.6","S.8")
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
colnames(dataTMM) = c("ID","S.1","S.2","S.3","S.4","S.5","S.6","S.8")
dataTMM = as.data.frame(dataTMM)
saveRDS(dataTMM, "dataTMM.rds")

# Make dataTMM2 object (see if get same TMM normalization factors)
y <- DGEList(counts=data, group= rep(1,7))
y <- calcNormFactors(y, method="TMM")
identical(tmmSizes, y[[2]]$norm.factors)
identical(data, as.data.frame(y[[1]]))

