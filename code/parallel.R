#test multicore

rm(list=ls())

library(parallel)

dats=list(n=10,k=5)

simdat = function(...){
  #include a pause to see if it works
  
res=list(n=10,k=5)
  
  #dats=c(...)
  dats=list(n=10,k=5)
  n=dats[['n']]
  k=dats[['k']]
  
  d=list(matrix(rnorm(n*k,mean=0,sd=1),n,k))
  
  
}

simdat2 = function(...){
  
  d=list(matrix(rnorm(n*k,mean=0,sd=1),n,k))
  
}

#xtest = simdat(n=5,k=10)
xtest = simdat()

cl <- makeCluster(mc <- getOption("cl.cores", 2))
cltst=do.call(c,parLapply(cl=cl,1:2,simdat))

#alt
#res=clusterEvalQ(cl,dats)
#cl2=do.call(c,res)
simdat2(cl2)

stopCluster(cl)