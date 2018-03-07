



for (i in 1:iter){
  #print(i)
  #simulate beta from mvn
  b[i,]=pars+t(rnorm(length(pars),mean=0,sd=1))%*%chol(s2[i]*xtxi)
  yhat[i,] = x %*% b[i,]
  sse = sum((y-(yhat[i,]))^2)
  sst = sum((y-mean(y))^2)
  r2[i] = 1-(sse/sst)
  rmse[i] = sqrt(sse/n)
  ll[i]=sum(dnorm(y,mean=yhat[i,],sd=s2[i],log=TRUE))
}

#use internal loop iterators which save some time

#prepare covariance
sdraw = lapply(s2,FUN=function(s) chol(s*xtxi))
#draw betas
bb = lapply(sdraw,FUN=function(sigma) pars + t(rnorm(length(pars),mean=0,sd=1))%*%sigma)
yhat = lapply(bb,FUN=function(beta) x %*% t(beta))
sse = lapply()

