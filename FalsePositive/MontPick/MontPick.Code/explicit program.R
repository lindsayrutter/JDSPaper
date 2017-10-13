source("input.pickrell.R")

packageRunning <- 
function(method, count.data, cond, n1=3, n2=3, nsim=100, nseed=200) {
  #record the time
  start.time <- proc.time()
  # separate data for simulation
  sub.data1 <- count.data[which(cond == unique(cond)[1])]
  sub.data2 <- count.data[which(cond == unique(cond)[2])]
  # choose method
  method <- method
  # simulation 
  source("sampling.r")
  get.name <- paste(method,"simul.R")
  source(get.name) 
  print(get.name)
  func.name <- paste0(method, ".simulation")
  tmp <- get(func.name)
  simulPtable1 <- tmp(sub.data1, n1, n2, nsim,nseed)
  print(proc.time()-start.time)
  simulPtable2 <- tmp(sub.data2, n1, n2, nsim,nseed)
  print(proc.time()-start.time)
  # norminal p table
  func.group <- paste0(method, ".group")
  temp <- get(func.group)
  normPtable <- temp(sub.data1, sub.data2, n1, n2, nsim,nseed)
  # calculate the corrected pvalue table 
  source("calcPower.r")
  if (method == "DESeq" || method == "DESeq2") {
  	  simulPtable1[is.na(simulPtable1)] <- 1
      simulPtable2[is.na(simulPtable2)] <- 1
      normPtable[is.na(normPtable)] <- 1
  }
  correctP <- calcPower(simulPtable1, simulPtable2, normPtable, nsim)
  #time recorded
  elapsed.time <- proc.time() - start.time
  print(elapsed.time)
  #save
  save(simulPtable1, simulPtable2, normPtable, correctP, file=paste(method," n=", n1, ".RData", sep = ''))
  print(paste("Your file has been saved in: ", method," n=", n1, ".RData", sep = ''))
  # return
  return(list(simPtable1=simulPtable1, simPtable2=simulPtable2, norminalP=normPtable, correctedP=correctP))
  
	
}




