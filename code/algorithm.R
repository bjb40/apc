###

###
#p. 563 Gelman
###

###
#basic linear model gibbs sampler (Scott)
###


lin_gibbs = function(y,x){
  iter = 1000
  
  r2=sse=sst=s2=matrix(1,iter)
  b= matrix(0,iter,ncol(x))
  yhat=matrix(0,length(y))
  xtxi = solve(t(x)%*%x)
  pars=coefficients(lm(y~x-1))
  
  #simulate sigma from inverse gamma marginal
  s2 = 1/rgamma(iter,nrow(x)-ncol(x)/2,.5*t(residuals(lm(y~x-1)))%*%residuals(lm(y~x-1)))
  
  #set ppd
  ppd = matrix(0,iter,length(y))
  
  #simulate beta from mvn
  #you can vetorize this!!
  for (i in 1:iter){
    b[i,]=pars+t(rnorm(length(pars),mean=0,sd=1))%*%chol(s2[i]*xtxi)
    yhat = x %*% b[i,]
    sse[i] = sum((y-yhat)^2)
    sst[i] = sum((y-mean(y))^2)
    r2[i] = 1-(sse[i]/sst[i])
  }
  
  ###BIC estimate for Bayes Factor (posterior probability weight)
  ###p. 135, Eq. 26 from Rafferty 1995 (SMR)
  n=length(y)
  bic=n*log(1-mean(r2))+(nrow(x)-1)*log(n)
  
  return(list(betas=b,sigma=s2,r2=r2,bic=bic))

}#end linear gibbs

#can add
res_gibbs=lin_gibbs(y=tdat$y2,x=model.matrix(~a+p,tdat))
res_lm=lm(y2~a+p,tdat)


##
tdat2=tst$newdat
gibbs2=lin_gibbs(y=tdat2$y2,model.matrix(~a+p+cohort,tdat2))


 