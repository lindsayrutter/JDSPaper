library(edgeR)
library(data.table)

data(soybean_cn)
data = soybean_cn
rownames(data) = data[,1]

y = DGEList(counts=data[,-1])
group = c(1,1,1,2,2,2,3,3,3)
y = DGEList(counts=y, group=group)
Group = factor(c(rep("S1",3), rep("S2",3), rep("S3",3)))
design <- model.matrix(~0+Group, data=y$samples)
colnames(design) <- levels(Group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

soybean_cn_metrics <- list()

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
    
    soybean_cn_metrics[[paste0(colnames(fit)[i], "_", colnames(fit)[j])]] <- lrt
  }
}

save(soybean_cn_metrics, file = "../data/soybean_cn_metrics.rda")
