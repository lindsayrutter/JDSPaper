library(limma)
library(edgeR)
print(paste("limma package version", packageVersion("limma")))
print(paste("edgeR package Version", packageVersion("edgeR")))

edgeR.robust.running <- function(count.data, cond) {
  group <- cond
  design <- model.matrix(~group, data=count.data)
  d <- DGEList(counts=count.data, group=group)
  d <- calcNormFactors(d)
  d <- estimateGLMRobustDisp(d, design)
  f <- glmFit(d, design=design)
  test <- glmLRT(f, coef=2)
  return(test)
}

edgeR.robust.simulation <- 
function(count.data, n1=3, n2=3, nsim=100, nseed) {
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim) {
    sub.samp <- sub.sample(count.data, n1, n2)
    out <- edgeR.robust.running(sub.samp, condition)
    store.p[, i] <- out$table$PValue
  }
  return(store.p)
}


edgeR.robust.group <- 
function(count.data1,count.data2, n1=3, n2=3, nsim=100, nseed) {
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim) {
    sub.samp <- sampling.groups(count.data1, count.data2, n1, n2)
    out <- edgeR.robust.running(sub.samp, condition)
    store.p[, i] <- out$table$PValue
  }
  return(store.p)
}
  
