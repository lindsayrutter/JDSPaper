source("input.pickrell.R")

numgenes <- dim(count.data)[1]
methods1 <- c("DESeq","DESeq2","edgeR","edgeR.robust","limma.voom")
#methods2 <- c("DESeq","DESeq2","edgeR","edgeR.robust","wald","lr")
methods2 <- c("DESeq","DESeq2","edgeR","edgeR.robust","lr")
methods3 <- c("limma.voom","limma.trans","ttest")



null.pv.plot <- function(n=3,numgenes,methods,fname){
  pdf(file=fname)
  pvec=c(.1,.01,.001,.0001,.00001,1e-6)
  num.methods <- length(methods)

  numsims <- 100
  numests <- num.methods
  numres <- length(pvec)
  numsig <- rep(0,numres)
  ratiosig <- rep(0,numres)
  fracsig <- rep(0,numres)

  ymin <- 1
  ymax <- 1e6
  expsig <- 2*pvec*numgenes*numsims
  plot(pvec,expsig,log="xy",ylim=c(ymin,ymax),xlab="P-Value",ylab="Number Significant",
       pch=0,lty=1,col="black",type="b",lwd=2)
  pchvec <- c(21:25,15)[1:num.methods]
  ltyvec <- (1:6)[1:num.methods]
  colvec <- c("red","blue","green","magenta","orange","brown")[1:num.methods]

  for (i in 1:num.methods){
  	mat.names <- paste(methods[i], " n=", n, ".RData", sep='')
  	df <- load(mat.names)
    pvs <- cbind(simulPtable1, simulPtable2)
    numsig <- rep(0,numres)
    for (j in 1:numres){
      numsig[j] <- sum(pvs < pvec[j],na.rm=T)
    }
print(numsig)

    lines(pvec,numsig,col=colvec[i],lty=ltyvec[i],lwd=2)
    points(pvec,numsig,col=colvec[i],pch=pchvec[i],lwd=2)    
  }
  ltext <- c("Expected Number Significant",methods)
  legend(3e-4,70,legend=ltext,col=c("black",colvec),
         pch=c(0,pchvec),lty=c(1:6,1),lwd=2)
  title(paste("Performance of Methods under the Null with n = ",n))
  dev.off()
}

power.plot <- function(n=3,numgenes,methods,corr=F,fname){
  pdf(file=fname)
  pvec=c(.1,.01,.001,.0001,.00001,1e-6)
  num.methods <- length(methods)

  numsims <- 100
  numests <- num.methods
  numres <- length(pvec)
  numsig <- rep(0,numres)
  ratiosig <- rep(0,numres)
  fracsig <- rep(0,numres)

  ymin <- 1
  ymax <- 1e6
  expsig <- pvec*numgenes*numsims
  plot(pvec,expsig,log="xy",ylim=c(ymin,ymax),xlab="P-Value",ylab="Number Significant",
       pch=0,lty=1,col="black",type="b",lwd=2)
  pchvec <- c(21:25,15)[1:num.methods]
  ltyvec <- (1:6)[1:num.methods]
  colvec <- c("red","blue","green","magenta","orange","brown")[1:num.methods]

  for (i in 1:num.methods){
    mat.names <- paste(methods[i], " n=", n, ".RData", sep='')
  	df <- load(mat.names)
  	pvs <- normPtable
        if (corr){
      	  pvs <- correctP
    	}
    numsig <- rep(0,numres)
    for (j in 1:numres){
      numsig[j] <- sum(pvs < pvec[j],na.rm=T)
    }
print(numsig)

    lines(pvec,numsig,col=colvec[i],lty=ltyvec[i],lwd=2)
    points(pvec,numsig,col=colvec[i],pch=pchvec[i],lwd=2)    
  }
  ltext <- c("Null Expected Number Significant",methods)
  legend(2e-4,70,legend=ltext,col=c("black",colvec),
         pch=c(0,pchvec),lty=c(1:6,1),lwd=2)
  title(paste("Power of Methods with n = ",n))
  dev.off()
}