require(MASS)

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

lr.running <- function(counts,cond,norm=T){
  if (norm) counts <- mednorm(counts)
  lr.func <- function(x,cond){
    var.pool <- sum(tapply(x,cond,var))
    if (var.pool < 1e-10){
      pv <- 1
    }
    else{
      pv <- tryCatch(anova(glm.nb(x~cond),glm.nb(x~1),test="Chisq")[2,8],error=etrap,silent=T)
    }
    return(pv)
  }
  pvs <- apply(counts,1,lr.func,cond=cond)
  return(data.frame(pval=pvs))
}

lr.simulation <- function(count.data, n1, n2, nsim, nseed,norm=T){
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim){
    sub.samp <- sub.sample(count.data, n1, n2)
    out <- lr.running(sub.samp, condition,norm=T)
    store.p[, i] <- out$pval
  }
  return(store.p)
}
    
lr.group <- 
function(sub.data1, sub.data2, n1=3, n2=3, nsim=100, nseed,norm=T){
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim){
    sub.samp <- sampling.groups(sub.data1, sub.data2, n1, n2)
    out <- lr.running(sub.samp, condition,norm=T)
    store.p[, i] <- out$pval
  }
  return(store.p)
}


