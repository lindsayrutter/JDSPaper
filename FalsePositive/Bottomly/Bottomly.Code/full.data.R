source("input.bottomly.R")
source("sampling.R")
n1 <- 10
n2 <- 11
n <- n1+n2
sub.data1 <- count.data[,1:n1]
sub.data2 <- count.data[,(n1+1):n]

source("DESeq simul.R")
DESeq.results <- DESeq.group(sub.data1, sub.data2, n1=n1, n2=n2, nsim=1, nseed=200)

source("DESeq2 simul.R")
DESeq2.results <- DESeq2.group(sub.data1, sub.data2, n1=n1, n2=n2, nsim=1, nseed=200)

source("edgeR simul.R")
edgeR.results <- edgeR.group(sub.data1, sub.data2, n1=n1, n2=n2, nsim=1, nseed=200)

source("edgeR.robust simul.R")
edgeR.robust.results <- edgeR.robust.group(sub.data1, sub.data2, n1=n1, n2=n2, nsim=1, nseed=200)

source("limma.voom simul.R")
limma.voom.results <- limma.voom.group(sub.data1, sub.data2, n1=n1, n2=n2, nsim=1, nseed=200)

all.results <- cbind(DESeq.results,DESeq2.results,edgeR.results,edgeR.robust.results,limma.voom.results)
colnames(all.results) <- c("DESeq","DESeq2","edgeR","edgeR.robust","limma.voom")
save(all.results,file="Bottomly.full.results.RData")
