library(limma)
library(edgeR)
print(paste("limma package version", packageVersion("limma")))
print(paste("edgeR package Version", packageVersion("edgeR")))

edgeR.running <- 
function(count.data, cond, mode.norm="TMM", dispersion.type="auto") {
  cds <- DGEList(count.data, group=cond)
  cds <- calcNormFactors(cds, method=mode.norm )
  cds <- estimateCommonDisp(cds)
  cds <- estimateTagwiseDisp(cds)
  et  <- exactTest(cds, dispersion=dispersion.type)
  return(et)
}

edgeR.simulation <- function(count.data, n1=3, n2=3, nsim=100, nseed) {
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim) {
    sub.samp <- sub.sample(count.data, n1, n2)
    out <- edgeR.running(sub.samp, condition)
    store.p[, i] <- out$table$PValue
  }
  return(store.p)
}
    

edgeR.group <- 
function(count.data1, count.data2, n1=3, n2=3, nsim=100, nseed) {
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim) {
    sub.samp <- sampling.groups(count.data1, count.data2, n1, n2)
    out  <- edgeR.running(sub.samp, condition)
    store.p[, i] <- out$table$PValue
  }
  return(store.p)
}
    