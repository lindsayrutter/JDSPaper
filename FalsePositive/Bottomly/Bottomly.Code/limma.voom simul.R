library(limma)
library(edgeR)
print(paste("limma package version", packageVersion("limma")))
print(paste("edgeR package Version", packageVersion("edgeR")))

limma.voom.running <- function(count.data, cond, method.norm="TMM" ) {
  dge <- DGEList(counts=count.data, genes=rownames(count.data), group=factor(cond))
  dge <- calcNormFactors(dge, method.norm)
  design <- model.matrix(~factor(cond))
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  return(fit)
}

limma.voom.simulation <- 
function(count.data, n1=3, n2=3, nsim=100, nseed) {
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim) {
    sub.samp <- sub.sample(count.data, n1, n2)
    out  <- limma.voom.running(sub.samp, condition)
    store.p[, i] <- out$p.value[, 2]
  }
  return(store.p)
}


limma.voom.group <- 
function(sub.data1, sub.data2, n1=3, n2=3, nsim=100, nseed) {
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim) {
    sub.samp <- sampling.groups(sub.data1, sub.data2, n1, n2)
    out  <- limma.voom.running(sub.samp, condition)
    store.p[, i] <- out$p.value[, 2]
  }
  return(store.p)
}