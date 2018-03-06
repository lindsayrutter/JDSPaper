source("plot.pickrell.R")
nvec <- c( 3, 4, 5, 10)
for (n in nvec){
#  fname <- paste("pickrell.null1.",n,".pdf",sep="")
#  null.pv.plot(n,numgenes=numgenes,methods=methods1,fname)
  fname <- paste("pickrell.null2.",n,".pdf",sep="")
  null.pv.plot(n,numgenes=numgenes,methods=methods2,fname)
#  fname <- paste("pickrell.null3.",n,".pdf",sep="")
#  null.pv.plot(n,numgenes=numgenes,methods=methods3,fname)
#  fname2 <- paste("pickrell.power3.",n,".pdf",sep="")
#  power.plot(n,numgenes=numgenes,methods=methods3,corr=F,fname=fname2)
}