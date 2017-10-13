library("DESeq2")
print(paste("DESeq2 package version", packageVersion("DESeq2")))

DESeq2.running <- 
function(count.data,cond, test.type="Wald", fit.type="parametric"){   
  cond <- factor(cond)
  dds <- DESeqDataSetFromMatrix(countData=count.data, design=~cond, colData=as.data.frame(cond))
  dds <- DESeq(dds, test=test.type, fitType=fit.type, quiet=T)
  res <- results(dds)
  return(res$pvalue)
}

DESeq2.simulation <- 
function(count.data, n1=3, n2=3, nsim=100, nseed) {
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim) {
    sub.samp <- sub.sample(count.data, n1, n2)
    out <- DESeq2.running(sub.samp, condition)
    store.p[, i] <- out
  }
  return(store.p)
}

    
DESeq2.group <- 
function(sub.data1, sub.data2, n1=3, n2=3, nsim=100, nseed) {
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim) {
    sub.samp <- sampling.groups(sub.data1, sub.data2, n1, n2)
    out <- DESeq2.running(sub.samp, condition)
    store.p[, i] <- out
  }
  return(store.p)
}



