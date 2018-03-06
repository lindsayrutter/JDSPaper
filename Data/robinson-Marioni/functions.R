plotReverseCumDist <- function(x, xlab="Tag count X", ylab="# tags >= X", add=FALSE, ...) {
  v <- ecdf(x)
  matplot( knots(v), (1-v(knots(v)))*sum(D[,1]), log="xy", xlab=xlab, ylab=ylab, add=add, ... )
}

generateDataset2 <- function(commonTags=15000, uniqueTags=c(1000,3000), group=c(1,2), libLimits=c(.9,1.1)*1e6, empiricalDist=NULL, lengthDist=NULL, pDifferential=.05, pUp=.5, foldDifference=2, nreps=c(2,2)) {
  
  # some checks
  stopifnot( length(group) == length(uniqueTags) )
  stopifnot( length(group) == length(nreps) )
  stopifnot( length(empiricalDist) == length(lengthDist) )
  group <- as.factor(rep(group,nreps))
  stopifnot( nlevels(group) == 2 ) 
  
  print(group)
  
  #exampleCounts <- empiricalDist/lengthDist
  exampleCounts <- empiricalDist
  exampleLambda <- exampleCounts/sum(exampleCounts)
  exampleIds <- seq_len( length(empiricalDist) )
  
  # set up libraries
  nLibraries <- sum( nreps )
  libSizes <- runif(nLibraries, min=libLimits[1], max=libLimits[2] )
  
  # vector of starts/stops for the unique Tags
  en <- commonTags + cumsum(uniqueTags)
  st <- c(commonTags+1,en[-nLibraries]+1)
  
  # create matrix of LAMBDA(=relative expression levels)
  LAMBDA <- matrix(0, nrow=max(en), ncol=nLibraries)
  
  ID <- rep(0, max(en))
  ID[1:commonTags] <- sample(exampleIds, commonTags, replace=TRUE)
  LAMBDA[1:commonTags,] <- exampleLambda[ ID[1:commonTags] ]
  
  # set unique tag totals
  for(i in 1:length(nreps))
    if(uniqueTags[i] > 0) {
      ID[st[i]:en[i]] <- sample(exampleIds, uniqueTags[i], replace=TRUE)
      LAMBDA[st[i]:en[i],group==levels(group)[i]] <- exampleLambda[ ID[st[i]:en[i]] ]
    }
  
  g <- group == levels(group)[1]
  ind <- seq_len(floor(pDifferential*commonTags))
  if(length(ind)>0) {
    fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
    LAMBDA[ind,g] <- LAMBDA[ind,g]*exp(log(foldDifference)/2*fcDir)
    LAMBDA[ind,!g] <- LAMBDA[ind,!g]*exp(log(foldDifference)/2*(-fcDir))
  }
  
  sampFactors <- colSums(LAMBDA)
  
  sampFactorsM <- outer( rep(1,max(en)), sampFactors )
  libSizesM <- outer(  rep(1,max(en)), libSizes )
  
  # create observed means
  MEAN <- LAMBDA / sampFactorsM * libSizesM  # to get the totals to sum to 1
  
  # sample observed data (column sums will be *close* to set library sizes)
  DATA <- matrix(0, nr=nrow(LAMBDA), ncol=nLibraries)
  DATA <- matrix(rpois(length(MEAN), lambda=MEAN),ncol=nLibraries)
  
  trueFactors <- colSums(MEAN[1:commonTags,])
  trueFactors <- trueFactors/trueFactors[1]
  
  colnames(DATA) <- paste(paste("group",group,sep=""),1:ncol(DATA),sep=".")
  
  list(DATA=DATA, LAMBDA=LAMBDA, MEAN=MEAN, trueFactors=trueFactors, group=group, libSizes=libSizes,  
       differentialInd=c(ind,(commonTags+1):nrow(DATA)), commonInd=1:commonTags, ID=ID, length=lengthDist[ID])
}

takeSubset <- function(obj, subsetInd) {
  allInd <- 1:nrow(obj$DATA)
  commonInd <- allInd %in% obj$commonInd
  differentialInd <- allInd %in% obj$differentialInd
  list(DATA=obj$DATA[subsetInd,], LAMBDA=obj$LAMBDA[subsetInd,], trueFactors=obj$trueFactors, group=obj$group, 
       libSizes=obj$libSizes, differentialInd=which(differentialInd[subsetInd]), commonInd=which(commonInd[subsetInd]),
       ID=obj$ID[subsetInd], length=obj$length[subsetInd])
}

generateDataset <- function(commonTags=15000, uniqueTags=c(1000,3000), group=c(1,2), libLimits=c(.9,1.1)*1e6, empiricalDist=NULL, randomRate=1/100, pDifferential=.05, pUp=.5, foldDifference=2) {
  
  # some checks
  group <- as.factor(group)
  stopifnot( length(group) == length(uniqueTags) )
  #stopifnot( length(group) == 2 ) # code below only works for 2 samples
  stopifnot( nlevels(group) == 2 ) 
  
  # define where to take random sample from (empirical distribution OR random exponential)
  if(is.null(empiricalDist))
    exampleCounts <- ceiling(rexp(commonTags,rate=randomRate))
  else
    exampleCounts <- empiricalDist
  
  exampleLambda <- exampleCounts/sum(exampleCounts)
  
  # set up libraries
  nLibraries <- length(uniqueTags)
  libSizes <- runif(nLibraries, min=libLimits[1], max=libLimits[2] )
  
  # vector of starts/stops for the unique Tags
  en <- commonTags + cumsum(uniqueTags)
  st <- c(commonTags+1,en[-nLibraries]+1)
  
  # create matrix of LAMBDA(=relative expression levels)
  LAMBDA <- matrix(0, nrow=max(en), ncol=nLibraries)
  LAMBDA[1:commonTags,] <- sample(exampleLambda, commonTags, replace=TRUE)
  
  # set unique tag totals
  for(i in 1:nLibraries)
    if(uniqueTags[i] > 0)
      LAMBDA[st[i]:en[i],i] <- sample(exampleLambda, uniqueTags[i])    
  
  ind <- seq_len(floor(pDifferential*commonTags))
  g <- group == levels(group)[1]
  if(length(ind)>0) {
    fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
    LAMBDA[ind,g] <- LAMBDA[ind,!g]*exp(log(foldDifference)/2*fcDir)
    LAMBDA[ind,!g] <- LAMBDA[ind,!g]*exp(log(foldDifference)/2*(-fcDir))
  }
  
  sampFactors <- colSums(LAMBDA)
  
  sampFactorsM <- outer( rep(1,max(en)), sampFactors )
  libSizesM <- outer(  rep(1,max(en)), libSizes )
  
  # create observed means
  MEAN <- LAMBDA / sampFactorsM * libSizesM
  
  # sample observed data (column sums will be *close* to set library sizes)
  DATA <- matrix(rpois(length(MEAN), lambda=MEAN),ncol=nLibraries)
  
  trueFactors <- colSums(MEAN[1:commonTags,])
  trueFactors <- trueFactors/trueFactors[1]
  list(DATA=DATA, LAMBDA=LAMBDA, MEAN=MEAN, trueFactors=trueFactors, group=group, libSizes=libSizes,  differentialInd=c(ind,(commonTags+1):nrow(DATA)), 
       commonInd=1:commonTags)
}

calcFactor <- function(obs, ref, trim=.45) {
  logR <- log2(obs/ref)
  fin <- is.finite(logR)
  2^mean(logR[fin],trim=trim)
}

Poisson.model <- function(MA,group1,group2){
  
  require(limma)
  Poisson.glm.pval <- vector()
  Fold.changes <- vector()
  
  CS <- colSums(MA$M[,c(group1,group2)])
  
  for (i in 1:(nrow(MA))){
    S1 <- MA$M[i,group1] 
    S2 <- MA$M[i,group2] 
    In <- c(S1,S2)
    sample.f <- factor(c(rep(1,length(group1)),rep(2,length(group2))))
    In <- as.vector(unlist(In))
    GLM.Poisson <- glm(In ~ 1 + sample.f + offset(log(CS)),family=poisson)
    Poisson.glm.pval[i] <- anova(GLM.Poisson,test="Chisq")[5][2,1]
    Fold.changes[i] <- exp(GLM.Poisson$coefficients[1])/(exp(GLM.Poisson$coefficients[1]+GLM.Poisson$coefficients[2]))
  }
  
  output <- matrix(ncol=2,nrow=nrow(MA$M))
  output[,1] <- Poisson.glm.pval
  output[,2] <- Fold.changes
  output <- as.data.frame(output)
  names(output) <- c("pval","FC")
  output
}

Poisson.model.new <- function(countMatrix,group1,group2, ref=1, calcFactor=TRUE){
  
  Poisson.glm.pval <- vector()
  Fold.changes <- vector()
  
  props <- countMatrix / outer( rep(1,nrow(countMatrix)), colSums(countMatrix) )
  
  refS <- colSums(countMatrix[,c(group1,group2)])
  
  if( calcFactor ) {
    require(edgeR)
    CS <- calcNormFactors(countMatrix[,c(group1,group2)])
  } else {
    CS <- rep(1,length(group1)+length(group2))
  }
  
  offsets <- log(CS)+log(refS)
  
  sample.f <- factor(c(rep(1,length(group1)),rep(2,length(group2))))
  
  for (i in 1:(nrow(countMatrix))){
    S1 <- countMatrix[i,group1] 
    S2 <- countMatrix[i,group2] 
    In <- c(S1,S2)
    In <- as.vector(unlist(In))
    GLM.Poisson <- glm(In ~ 1 + sample.f + offset(offsets),family=poisson)
    Poisson.glm.pval[i] <- anova(GLM.Poisson,test="Chisq")[5][2,1]
    Fold.changes[i] <- exp(GLM.Poisson$coefficients[1])/(exp(GLM.Poisson$coefficients[1]+GLM.Poisson$coefficients[2]))
    if(i %% 100==0) cat(".")
  }
  cat("\n")
  
  #output <- matrix(ncol=2,nrow=nrow(countMatrix))
  #output[,1] <- Poisson.glm.pval
  #output[,2] <- Fold.changes
  #output <- as.data.frame(output)
  #names(output) <- c("pval","FC")
  
  list(stats=data.frame(pval=Poisson.glm.pval, FC=Fold.changes),offsets=offsets,factors=CS)
}


exactTestPoisson <- function(dataMatrix, meanMatrix, group1Ind, group2Ind, verbose=TRUE) {
  
  y1 <- rowSums(dataMatrix[,group1Ind])
  y2 <- rowSums(dataMatrix[,group2Ind])
  m1 <- rowSums(meanMatrix[,group1Ind])
  m2 <- rowSums(meanMatrix[,group2Ind])
  
  N <- rowSums( dataMatrix[,c(group1Ind,group2Ind)] )
  
  pvals <- rep(NA, nrow(dataMatrix))
  
  for (i in 1:length(pvals)) {
    v <- 0:N[i]
    p.top <- dpois(v, lambda=m1[i]) * dpois(N[i]-v, lambda=m2[i])
    p.obs <- dpois(y1[i], lambda=m1[i]) * dpois(y2[i], lambda=m2[i])
    p.bot <- dpois(N[i], lambda=m1[i]+m2[i])
    keep <- p.top <= p.obs
    pvals[i] <- sum(p.top[keep]/p.bot)
    if (verbose)
      if (i%%1000 == 0)
        cat(".")
  }
  if (verbose)
    cat("\n")
  
  pvals
  
}

calcFactorRLM <- function(obs, ref, logratioTrim=.20, sumTrim=0.01) {
  
  if( all(obs==ref) )
    return(1)
  
  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))         # log ratio of expression, accounting for library size
  p0 <- obs/nO
  pR <- ref/nR
  
  x <- log2(p0)-log2(pR)
  x <- x[ !is.na(x) & is.finite(x) ]
  
  r <- rlm(x~1, method="MM")
  2^r$coef
  
}

calcFactorWeighted <- function(obs, ref, logratioTrim=.3, sumTrim=0.05) {
  
  if( all(obs==ref) )
    return(1)
  
  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))         # log ratio of expression, accounting for library size
  absE <- log2(obs/nO) + log2(ref/nR)     # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref  # estimated asymptotic variance
  
  fin <- is.finite(logR) & is.finite(absE)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  
  # taken from the original mean() function
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  
  keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  2^( sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE) )
  
}

calcFactor2 <- function(obs, ref) {
  logR <- log2(obs/ref)
  fin <- is.finite(logR)
  d<-density(logR,na.rm=TRUE)
  2^d$x[which.max(d$y)]
}


fdPlot <- function( score, indDiff, add=FALSE, xlab="Number of Genes Selected", 
                    ylab="Number of False Discoveries", lwd=4, type="l", ... ) {
  o <- order(score)
  w <- o %in% indDiff
  x <- 1:length(indDiff)
  y <- cumsum(!w[indDiff])
  matplot(x, y, xlab=xlab, ylab=ylab, lwd=lwd, type=type, add=add, ... )
}