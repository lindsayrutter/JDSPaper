mednorm <- function(mat){
  cmeds <- apply(mat,2,median)
  medmed <- median(cmeds)
  appnorm <- function(x){
    return(medmed*x/cmeds)
  }
  norm.mat <- t(apply(mat,1,appnorm))
  return(norm.mat)
}

etrap <- function(econd){ 
  return(1)
}

ttest.running <- function(counts,cond,norm=T,var.equal=F){
  if (norm) counts <- mednorm(counts)
  t.func <- function(x,cond){
    var.pool <- sum(tapply(x,cond,var))
    if (var.pool < 1e-10){
      pv <- 1
    }
    else{
      pv <- tryCatch(t.test(x~cond,var.equal=var.equal)$p.value,error=etrap,silent=T)
    }
  }
  pvs <- apply(counts,1,t.func,cond=cond)
  return(data.frame(pval=pvs))
}


ttest.simulation <- function(count.data, n1, n2, nsim, nseed,norm=T,var.equal=F){
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim){
    sub.samp <- sub.sample(count.data, n1, n2)
    out <- ttest.running(sub.samp, condition,norm=T,var.equal=F)
    store.p[, i] <- out$pval
  }
  return(store.p)
}
    
ttest.group <- 
function(sub.data1, sub.data2, n1=3, n2=3, nsim=100, nseed,norm=T,var.equal=F){
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim){
    sub.samp <- sampling.groups(sub.data1, sub.data2, n1, n2)
    out <- ttest.running(sub.samp, condition,norm=T,var.equal=F)
    store.p[, i] <- out$pval
  }
  return(store.p)
}


