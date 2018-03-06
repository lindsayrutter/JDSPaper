# This is the data from reference 18 of TMM Robinson paper.
# It is the data that has been normalized by library size, created in makeDataLibTMM.R

data = readRDS("dataLib.rds")

baseOutDir = "/Users/lindz/JDSPaper/Data/robinson-Cloonan/libNorm"

# Obtain ESEB values
outDir = paste0(baseOutDir, "/ESEB")
dataSel <- data
dataSel[,c(2:ncol(dataSel))] = log(dataSel[,c(2:ncol(dataSel))]+1)
plotScatterStatic(dataSel, outDir = outDir)

