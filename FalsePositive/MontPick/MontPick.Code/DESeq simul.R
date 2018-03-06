library(DESeq)
print(paste("DESeq Package Version", packageVersion("DESeq")))

DESeq.running <- 
function(count.data, cond, method="pooled", sharingmode="maximum"){

  ct <- newCountDataSet(count.data, cond)
  ct <- estimateSizeFactors(ct)
  ct <- estimateDispersions(ct, method, sharingmode)
  resi <- nbinomTest(ct, as.character(unique(cond)[1]), as.character(unique(cond)[2]))
  return(resi)
}

DESeq.simulation <- function(count.data, n1, n2, nsim, nseed){
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim){
    sub.samp <- sub.sample(count.data, n1, n2)
    out <- DESeq.running(sub.samp, condition)
    store.p[, i] <- out$pval
  }
  return(store.p)
}
    
DESeq.group <- 
function(sub.data1, sub.data2, n1=3, n2=3, nsim=100, nseed){
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim){
    sub.samp <- sampling.groups(sub.data1, sub.data2, n1, n2)
    out <- DESeq.running(sub.samp, condition)
    store.p[, i] <- out$pval
  }
  return(store.p)
}


