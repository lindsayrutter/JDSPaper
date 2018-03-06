# This is the data from Section 4.5 of edgeR vignette
# https://mayoclinic.pure.elsevier.com/en/publications/tumor-transcriptome-sequencing-reveals-allelic-expression-imbalan

library(tweeDEseqCountData)
library(data.table)
data(pickrell1)
Counts <- exprs(pickrell1.eset)
data <- as.data.frame(Counts)
Gender <- pickrell1.eset$gender

# Obtain 3 male 3 female
data <- data[,c(1:5,7)]
setDT(data, keep.rownames = TRUE)[]
colnames(data) <- c("ID","m.1","m.2","f.1","m.3","f.2","f.3")
data <- as.data.frame(data)
data = data[,c(1,2,3,5,4,6,7)]
data[,c(2:ncol(data))] = log(data[,c(2:ncol(data))]+1)
plotScatterStatic(data)



