data <- read.table("../Bottomly.Data/bottomly_count_table.txt",header=T)
count.data <- data[, -1]
rownames(count.data) <- data[, 1]
min.count <- 21
count.data <- count.data[which(rowSums(count.data) >= min.count), ]
data2 <- read.table("../Bottomly.Data/bottomly_phenodata.txt", header=T)
condition <- factor(as.character(data2$strain))
