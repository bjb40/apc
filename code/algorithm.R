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
  
  ll=r2=s2=matrix(1,iter)
  b= matrix(0,iter,ncol(x))
  yhat=matrix(0,length(y))
  xtxi = solve(t(x)%*%x)
  pars=coefficients(lm(y~x-1))
  
  #simulate sigma from inverse gamma marginal
  s2 = 1/rgamma(iter,nrow(x)-ncol(x)/2,.5*t(residuals(lm(y~x-1)))%*%residuals(lm(y~x-1)))
  
  #set ppd
  ppd = matrix(0,iter,length(y))
  
  #simulate beta from mvn
  #if you can vetorize this you will speed up a lot !!
  for (i in 1:iter){
    b[i,]=pars+t(rnorm(length(pars),mean=0,sd=1))%*%chol(s2[i]*xtxi)
    sse = sum((y-(x %*% b[i,]))^2)
    sst = sum((y-mean(y))^2)
    r2[i] = 1-(sse/sst)
    ll[i]=sum(dnorm(y,mean=y-(x%*% b[i,]),sd=s2[i],log=TRUE))
  }
  
  colnames(b) = colnames(x)
  ###BIC estimate for Bayes Factor (posterior probability weight)
  ###p. 135, Eq. 26 from Rafferty 1995 (SMR)
  n=length(y)
  bic_prime=n*log(1-mean(r2))+(ncol(x)-1)*log(n)
  bic=2*mean(ll)+ncol(x)*log(n)
  
  return(list(betas=b,sigma=s2,r2=r2,bic=bic,bic_prime=bic_prime,ll=ll))

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

for(age_w in 0:5){
  
  #reset dataframe
  x=tdat[,c('a','p','c')]
  
  if(age_w==0){
    x = x[!colnames(x) == 'a']
  } else{
    x$a=window(tdat$a,winlength=age_w)
  }

  for(period_w in 0:5){

    if(period_w==0){
      x = x[!colnames(x) == 'p']
    } else{
      x$p=window(tdat$p,winlength=period_w)
    }
    
    for(cohort_w in 0:5){

      if(cohort_w==0){
        x = x[!colnames(x) == 'c']
      } else{
        x$c=window(tdat$c,winlength=cohort_w)
      }
      #skip unidentified models
      comb_w = c(age_w,period_w,cohort_w) 
      if(all(comb_w ==1) | all(comb_w ==0)){next} 
 
      cat('\n\nEstimating model',length(allmods),
            'with windows:',
            '\n\tage    ',age_w,
            '\n\tperiod ',period_w,
            '\n\tcohort ',cohort_w,
            '\n\n')
            
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


##post-processing --- best model
bics=unlist(lapply(effects,FUN=function(x) x$bic))
bics_prime=unlist(lapply(allmods,FUN=function(x) x$bic_prime))
r2=unlist(lapply(allmods,FUN=function(x) mean(x$r2)))

best=which(bics==min(bics))

print(
  t(
    apply(allmods[[best]]$betas,2,FUN=function(x) c(eff=mean(x),sd=sd(x)))
  )
)

print(mean(allmods[[best]]$r2))

#cohort plots
#NOTE that the zero  may be a messed up order in the name!!

library(ggplot2)

load(paste0(datdir,'luo_sim_fits.RData'))

preds=as.data.frame(apply(effects[[best]]$c,2,mean)); colnames(preds)='est'
rng=apply(effects[[best]]$c,2,quantile,c(0.025,0.975))
preds$up = rng[2,]
preds$down = rng[1,]
preds$actual=pltdat$s1[order(pltdat$cohort)]

g=ggplot(preds,
         aes(y=est,x=1:nrow(preds))) + 
  geom_point() + 
  geom_errorbar(aes(ymax=up,ymin=down)) +
  geom_line(aes(y=actual))

print(g)

##post-processing -- model averaging
#averaging algorithm from ... eq 35 from Rafferty SMR

#w = exp(-.5*bics)/sum(exp(-.5*bics))
#w_prime=exp(-.5*bics_prime)/sum(exp(-.5*bics_prime))

#there are severe underflow problems.... since this is just a proportionality issue
#/changing proportion
prop=0.000005
w = exp(prop*bics)/sum(exp(prop*bics))

for(m in seq_along(effects)){
  effects[[m]]$w=w[m]
}

#weighted mean

for(l in seq_along(effects)){
  for(t in c('a','p','c')){
    if(is.null(effects[[l]][[t]])){
      effects[[l]][[t]] = 0
    }
  }
}

wtmn=function(l,var,w){
  #makes weighted mean
  #l is list results
  #var is charactername
  r=NA
  if(length(l[[var]])==1){r=0} else{
    r=apply(l[[var]],2,mean)*w
  }
  return(r)
}

preds$m_est =  rowSums(
  as.data.frame(
      lapply(effects,FUN=function(e)
            wtmn(e,'c',e$w))
))


wtquant=function(l,var,w,q){
  #takes weighted mean of quantile
  #l is list results
  #var is charactername
  #q is quantile
  r=NA
  if(length(l[[var]])==1){r=0} else{
    r=apply(l[[var]],2,quantile,q)*w
  }
  return(r)
}

preds$m_down=rowSums(
  as.data.frame(lapply(
  effects,FUN=function(e)
    wtquant(e,'c',e$w,0.025))
          )
)

preds$m_up=rowSums(
  as.data.frame(lapply(
  effects,FUN=function(e)
    wtquant(e,'c',e$w,0.975))
  )
)

preds$cohort=1:nrow(preds)
#predsm=preds[,! colnames(preds) %in% c('up','down')]
#pp=melt(predsm,id='cohort')

m=ggplot(preds,
         aes(y=m_est,x=cohort)) + 
         geom_point() + 
         geom_errorbar(aes(ymax=m_up,ymin=m_down)) +
         geom_line(aes(y=actual))

       


print(m)  