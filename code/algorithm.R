###

#clear cache
rm(list=ls())
source('config~.R')


#load test data
load(paste0(datdir,'testdat.RData'))

###
#p. 563 Gelman
###

###
#basic linear model gibbs sampler (Scott)
###


lin_gibbs = function(y,x){
  iter = 1000
  
  r2=s2=matrix(1,iter)
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
    sse = sum((y-(x %*% b[i,]))^2)
    sst = sum((y-mean(y))^2)
    r2[i] = 1-(sse/sst)
  }
  
  colnames(b) = colnames(x)
  ###BIC estimate for Bayes Factor (posterior probability weight)
  ###p. 135, Eq. 26 from Rafferty 1995 (SMR)
  n=length(y)
  bic=n*log(1-mean(r2))+(nrow(x)-1)*log(n)
  
  return(list(betas=b,sigma=s2,r2=r2,bic=bic))

}#end linear gibbs


#@@@@@@@@@@@
#sampling from model subspace
#@@@@@@@@@@@

####
#nested for-loops; windows of 1-5 spaces each (including non-inclusion... (i.e. 0))
####

#calculate time; based on 1 second per calculation
models=6^3
cat('Estimated hours:',models/60/60)

y=tdat$y1
tdat$c=tdat$p-tdat$a

allmods=list() #may run into size constraints/may need to limit to best mods... 
effects=list()
tm=Sys.time()
avtm=0

for(age_w in 0:6){
  
  #reset dataframe
  x=tdat[,c('a','p','c')]
  
  if(age_w==0){
    x = x[!colnames(x) == 'a']
  } else{
    x$a=window(tdat$a,winlength=age_w)
  }

  for(period_w in 0:6){

    if(period_w==0){
      x = x[!colnames(x) == 'p']
    } else{
      x$p=window(tdat$p,winlength=period_w)
    }
    
    for(cohort_w in 0:6){

      if(cohort_w==0){
        x = x[!colnames(x) == 'c']
      } else{
        x$c=window(tdat$c,winlength=cohort_w)
      }
      
      if(age_w %in% c(0,1) & period_w %in% c(0,1)){next} ##this leaves out a bunch of models
 
      cat('\n\nEstimating model',length(allmods),'\n\n')
            
      #add esitmate
      mnum = length(allmods)+1
      m = allmods[[mnum]] = lin_gibbs(y=y,x=model.matrix(~.,x))
      
      #need post processing (b_hat) here through scopedummy
      #consider limiting base on occam's window...
      blockdat=lapply(x,scopedummy)
      predat=lapply(blockdat,FUN=function(x) 
        model.matrix(~.,data=as.data.frame(x)))

      effects[[mnum]] = list()
      
      betas=list()
      for(eff in names(predat)){      
        betas[[eff]] = m$betas[,grepl(paste0('Intercept|',eff),colnames(m$betas))]
        effects[[mnum]][[eff]] = t(predat[[eff]] %*% t(betas[[eff]]))
        colnames(effects[[mnum]][[eff]]) = paste0(eff,unique(tdat[,eff]))
      }
      
      effects[[mnum]]$bic=m$bic

#      if(length(allmods)%%10==0){
        avtm=(avtm*(length(allmods)-1)+Sys.time()-tm)/length(allmods)
        cat('Average model time:',avtm,'\n\n')
        tm=Sys.time()
#      }
    }#end cohort loop
  
  }#end period loop   
} #end age loop


#post-processing --- model-averaging based on relative BICs...


 