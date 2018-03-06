library(edgeR)
library(data.table)

data(soybean_cn)
data = soybean_cn
rownames(data) = data[,1]

# Switch and relabel
data <- data[,c(1,2,3,8,5,6,7,4,9,10)]
colnames(data) <- c("ID","S1.1","S1.2","S1.3","S2.1","S2.2","S2.3","S3.1","S3.2", "S3.3")

y = DGEList(counts=data[,-1])
group = c(1,1,1,2,2,2,3,3,3)
y = DGEList(counts=y, group=group)
Group = factor(c(rep("S1",3), rep("S2",3), rep("S3",3)))
design <- model.matrix(~0+Group, data=y$samples)
colnames(design) <- levels(Group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

soybean_cn_switch13_metrics <- list()

for (i in 1:(ncol(fit)-1)){
  for (j in (i+1):ncol(fit)){
    contrast=rep(0,ncol(fit))
    contrast[i]=1
    contrast[j]=-1
    lrt <- glmLRT(fit, contrast=contrast)
    lrt <- topTags(lrt, n = nrow(y[[1]]))[[1]]
    
    setDT(lrt, keep.rownames = TRUE)[]
    colnames(lrt)[1] = "ID"
    lrt <- as.data.frame(lrt)
    
    soybean_cn_switch13_metrics[[paste0(colnames(fit)[i], "_", colnames(fit)[j])]] <- lrt
  }
}

save(soybean_cn_switch13_metrics, file = "soybean_cn_switch13_metrics.rda")
