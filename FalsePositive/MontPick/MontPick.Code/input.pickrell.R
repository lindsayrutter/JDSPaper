data <- read.table("montpick_count_table.txt",header=T)
count.data <- data[, -1]
rownames(count.data) <- data[, 1]
min.count <- 129
count.data <- count.data[which(rowSums(count.data) >= min.count), ]
data2 <- read.table("montpick_phenodata.txt", header=T)
condition <- factor(as.character(data2$population))
