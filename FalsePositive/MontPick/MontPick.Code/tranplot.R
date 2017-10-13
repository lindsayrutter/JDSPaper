require(limma)
#source("input.pickrell.R")

mv.calc <- function(x,tmat,cond){
  require(matrixStats)
  levs <- unique(cond)
  if (x < -2) tmat <- as.matrix(tmat)
  if (x > -2 & x < 0) tmat <- as.matrix(log(tmat))
  if (x > 0) tmat <- as.matrix(tran1(tmat,x))
  means1 <- rowMeans(tmat[,cond==levs[1]])
  means2 <- rowMeans(tmat[,cond==levs[2]])
  vars1 <- rowVars(tmat[,cond==levs[1]])
  vars2 <- rowVars(tmat[,cond==levs[2]])
  vars <- (vars1*59+vars2*68)/127
  mv <- mean(vars)
  vv <- var(vars)
  nu <- 127
  evv <- 2*mv^2/nu
  ahat <- (mv^2*(1-2/nu)+2*vv)/(vv-2*mv^2/nu)
  bhat <- 1/(mv*(ahat-1))
  print(c(ahat,bhat))
  print(1/(ahat*bhat))
  print(mean(1/vars))
  return(c(mv,vv,evv))
}


tran1 <- function(y,c1)
{
  ty <- log(y+sqrt(y^2+y/c1^2)+1/(2*c1^2))
  return(ty)
}

mednorm <- function(mat){
  cmeds <- apply(mat,2,median)
  medmed <- median(cmeds)
  appnorm <- function(x){
    return(medmed*x/cmeds)
  }
  norm.mat <- t(apply(mat,1,appnorm))
  return(norm.mat)
}

tran.plot <- function(x,tmat,cond,title1){
  require(matrixStats)
  levs <- unique(cond)
  if (x < -2) tmat <- as.matrix(tmat)
  if (x > -2 & x < 0) tmat <- as.matrix(log(tmat))
  if (x > 0) tmat <- as.matrix(tran1(tmat,x))
  means1 <- rowMeans(tmat[,cond==levs[1]])
  means2 <- rowMeans(tmat[,cond==levs[2]])
  vars1 <- rowVars(tmat[,cond==levs[1]])
  vars2 <- rowVars(tmat[,cond==levs[2]])
  means <- c(means1,means2)
  vars <- c(vars1,vars2)
  mcut <- tran1(0.2,x)
  msub <- means > mcut
  means <- means[msub]
  vars <- vars[msub]
  if (x < -2){
    plot(means,vars,log="xy",xlab="Mean",ylab="Variance",main=title1)
#    lines(lowess(means,vars,f=.2),col="red",lwd=3)
#    abline(coef(lm(log(vars)~log(means))),col="red",lwd=3)
  }
  if (x > -2){ 
    plot(means,vars,xlab="Mean",ylab="Variance",main=title1)
    lines(lowess(means,vars,f=.2),col="red",lwd=3)
#    abline(coef(lm(vars~means)),col="red",lwd=3)
  }
  slope <- (coef(lm(vars~means))[2])
#  print(slope)
  return(slope)
}

filter1 <- function(count.data,cond,mlim){
  require(matrixStats)
  levs <- unique(cond)
  tmat <- as.matrix(count.data)
  means1 <- rowMeans(tmat[,cond==levs[1]])
  means2 <- rowMeans(tmat[,cond==levs[2]])
  msub <- means1 > mlim & means2 > mlim
  return(count.data[msub,])
}

tran.est <- function(count.data,cond){
  tmp <- uniroot(tran.reg,interval=c(0.1,1),tmat=as.matrix(count.data),cond=cond)
  return(tmp$root)
}

tran.reg <- function(x,tmat,cond){
#  start.time <- proc.time()
  require(matrixStats)
  levs <- unique(cond)
  tmat <- as.matrix(tran1(tmat,x))
  means1 <- rowMeans(tmat[,cond==levs[1]])
  means2 <- rowMeans(tmat[,cond==levs[2]])
  vars1 <- rowVars(tmat[,cond==levs[1]])
  vars2 <- rowVars(tmat[,cond==levs[2]])
  means <- c(means1,means2)
  vars <- c(vars1,vars2)
#  plot(means,vars)
#  print(x)
#  print(proc.time()-start.time)
  slope <- (coef(lm(vars~means))[2])
#  print(slope)
  return(slope)
}



limma.trans.running <- function(count.data, cond, method.norm="TMM" ) {
  c1 <- tran.est(count.data,cond)
  count.data <- mednorm(tran1(count.data,c1))
#  dge <- DGEList(counts=count.data, genes=rownames(count.data), group=factor(cond))
#  dge <- calcNormFactors(dge, method.norm)
  design <- model.matrix(~factor(cond))
#  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(count.data, design)
  fit <- eBayes(fit)
  return(fit)
}

limma.trans.simulation <- 
function(count.data, n1=3, n2=3, nsim=100, nseed) {
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim) {
    sub.samp <- sub.sample(count.data, n1, n2)
    out  <- limma.trans.running(sub.samp, condition)
    store.p[, i] <- out$p.value[, 2]
  }
  return(store.p)
}


limma.trans.group <- 
function(sub.data1, sub.data2, n1=3, n2=3, nsim=100, nseed) {
  set.seed(nseed)
  ngene <- nrow(count.data)
  store.p <- matrix(0, nrow=ngene, ncol=nsim)
  condition <- c(rep(0, n1), rep(1, n2))
  for (i in 1:nsim) {
    sub.samp <- sampling.groups(sub.data1, sub.data2, n1, n2)
    out  <- limma.trans.running(sub.samp, condition)
    store.p[, i] <- out$p.value[, 2]
  }
  return(store.p)
}