source("plot.bottomly.R")
nvec <- c( 3, 4, 5)
for (n in nvec){
  fname <- paste("bottomly.null1.",n,".pdf",sep="")
  null.pv.plot(n,numgenes=numgenes,methods=methods1,fname)
  fname <- paste("bottomly.null2.",n,".pdf",sep="")
  null.pv.plot(n,numgenes=numgenes,methods=methods2,fname)
  fname <- paste("bottomly.null3.",n,".pdf",sep="")
  null.pv.plot(n,numgenes=numgenes,methods=methods3,fname)
  fname2 <- paste("bottomly.power3.",n,".pdf",sep="")
  power.plot(n,numgenes=numgenes,methods=methods3,corr=F,fname=fname2)
}
