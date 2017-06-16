###

#clear cache
rm(list=ls())
source('config~.R')

#prelim
dv='y1'
actual='s1'

#load test data
#load(paste0(datdir,'testdat.RData'))
load(paste0(datdir,'nsim.RData'))

#@@@@@@@@@@@
#sampling from model subspace
#@@@@@@@@@@@

y=tdat[,dv]
tdat$c=tdat$p-tdat$a

allmods=list() #may run into size constraints/may need to limit to best mods... 
effects=xhats=ppd=list()
tm=Sys.time()
avtm=0

window.sample=function(var){
  #input is continuous of a,p,c
  #output is factor with uniform, random window constraints
  #see dirichelet, e.g. for alternative algorithms
  vals=unique(var)
  
  mean.wins=round(runif(1)*length(vals))
  winprob=1-(mean.wins/length(vals))
  partition=runif(length(vals)-1)>winprob
  
  #assign windows of 1 if all false
  if(!any(partition)){
    wins=window(var,winlength=1)
  } else{
    breaks=vals[which(partition==TRUE)]
    breaks=unique(c(min(var)-1,breaks,max(var)))
    wins=window(var,breaks=breaks)
  }
  
  return(wins)
}

#set of numbers of random samples
n.samples=25

#holder df for model summary data 
win = data.frame(a=numeric(), p=numeric(), c=numeric())
breaks=list(a=list(),p=list(),c=list())

for(s in 1:n.samples){

  #reset dataframe
  x=tdat[,c('a','p','c')]
  
  #draw random window samples
  x$a=window.sample(x$a)
  x$p=window.sample(x$p)
  x$c=window.sample(x$c)
  
  #collect model data
  nr=data.frame(a=length(levels(x$a)),
                p=length(levels(x$p)),
                c=length(levels(x$c)))
  win=rbind(win,nr)
  
  #collect breaks data
  breaks$a[[s]]=attr(x$a,'breaks')
  breaks$p[[s]]=attr(x$p,'breaks')
  breaks$c[[s]]=attr(x$c,'breaks')
  

      cat('\n\nEstimating model',length(allmods),
            'with windows:',
            '\n\tage    ',nr$a,
            '\n\tperiod ',nr$p,
            '\n\tcohort ',nr$c,
            '\n\n')
            
      #add esitmate
      mnum = length(allmods)+1
      m = allmods[[mnum]] = lin_gibbs(y=y,x=model.matrix(~.,x))
      
      #consider limiting base on occam's window...
      #how can I incorporate grandmeans into calculation of beta-hat?
      #also, some don't have this .... 
      grand.means=t(as.matrix(colSums(model.matrix(~.,x))/nrow(x)))
      
      blockdat=lapply(x,scopedummy)
      predat=lapply(blockdat,FUN=function(x) 
        model.matrix(~.,data=as.data.frame(x)))

      effects[[mnum]] = xhats[[mnum]] = list()
      
      betas=list()
      for(eff in names(predat)){
        #fix colnames
        colnames(predat[[eff]]) = sub('x',eff,colnames(predat[[eff]]))
        #input means for xhat
        
        xhat=grand.means[rep(seq(nrow(grand.means)), nrow(predat[[eff]])),]
        #replace specific with predat
        xhat[,colnames(predat[[eff]])] = predat[[eff]]
        xhats[[mnum]][[eff]] = xhat

        effects[[mnum]][[eff]] = t(xhat %*% t(m$betas))
        colnames(effects[[mnum]][[eff]]) = paste0(eff,unique(tdat[,eff]))
      }
      
      effects[[mnum]]$bic=m$bic
      
      #yhat
      #ppd[[mnum]] = xhat %*% t(as.matrix(allmods[[mnum]]$betas))
      #ytilde
      #ppd[[mnum]] = 
      
#      if(length(allmods)%%10==0){
        avtm=(avtm*(length(allmods)-1)+Sys.time()-tm)/length(allmods)
        cat('Average model time:',avtm,'\n\n')
        tm=Sys.time()

}#end sampling loop


##post-processing --- best model
bics=unlist(lapply(effects,FUN=function(x) x$bic))
bics_prime=unlist(lapply(allmods,FUN=function(x) x$bic_prime))
r2=unlist(lapply(allmods,FUN=function(x) mean(x$r2)))

best=which(bics_prime==min(bics_prime))

print(
  t(
    apply(allmods[[best]]$betas,2,FUN=function(x) c(eff=mean(x),sd=sd(x)))
  )
)

print(mean(allmods[[best]]$r2))

#cohort plots
#NOTE that the zero  may be a messed up order in the name!!

library(ggplot2)
library(gridExtra)

load(paste0(datdir,'luo_sim_fits.RData'))

best.plt = list()
preds = list()


#for(d in c('a','p','c')){
for(d in c('a','p','c')){
  if(effects[[best]][[d]]==0){
    effects[[best]][[d]]=matrix(0,1000,20)
    if(d == 'c'){
      effects[[best]][[d]]=cbind(effects[[best]][[d]],matrix(0,1000,19))
    }
  }
    
  preds[[d]]=as.data.frame(apply(effects[[best]][[d]],2,mean))
  colnames(preds[[d]])='est'
  rng=apply(effects[[best]][[d]],2,quantile,c(0.025,0.975))
  preds[[d]]$up = rng[2,]
  preds[[d]]$down = rng[1,]
  #s1-s4 are for scenarios --- needs to match with y1-y4
  #preds[[d]]$actual=pltdat[[d]]$s1[order(pltdat[[d]]$id)]
  pltdat = tdat %>%
    group_by_(d) %>% 
    summarize(s1=mean(s1)) %>%
    rename_(id=d)
  preds[[d]]$actual=pltdat[,'s1']
  preds[[d]]$id=pltdat$id
  
  
  #need to change axes and limits...actually better to melt into df 
  #and not to use grid.arrange... (same for below)
  calcplt=ggplot(preds[[d]],
    aes_string(y='est',x='id') ) + 
    geom_point() + 
    geom_errorbar(aes(ymax=up,ymin=down)) +
    geom_line(aes(y=actual))

  #print(calcplt)
  #Sys.sleep(2)
  
  best.plt[[d]]=calcplt
  
}

pdf(file=paste0(imdir,'best-fit.pdf'))

  pltrange=range(unlist(
    lapply(preds,function(x) range(x %>% select(-id)))
  ))
  
  grid.arrange(best.plt[['a']] + ylim(pltrange),
              best.plt[['p']] + ylim(pltrange),
              best.plt[['c']] + ylim(pltrange),
              ncol=3)

dev.off()

##post-processing -- model averaging
#averaging algorithm from ... eq 35 from Rafferty SMR
#http://scicomp.stackexchange.com/questions/1122/how-to-add-large-exponential-terms-reliably-without-overflow-errors
#http://stats.stackexchange.com/questions/249888/use-bic-or-aic-as-approximation-for-bayesian-model-averaging
#fixing overflow issue from suggestion above after logging to take difference

#w = exp(-.5*bics)/sum(exp(-.5*bics)))
#w_prime=exp(-.5*bics_prime)/sum(exp(-.5*bics_prime))

#sbics = bics[order(bics)][1:10]

k=min(bics)
d=-.5*(bics-k)
w=exp(d)/sum(exp(d))

k = min(bics_prime)
d=-.5*(bics_prime-k)
w_prime=exp(d)/sum(exp(d))


#add weight to apc windows dataframe
win$wt=w; win$w_prime=w_prime; win$r2=r2; win$bic=bics; win$bic_prime=bics_prime
win$rmse = unlist(lapply(allmods,FUN=function(x) mean(x$rmse))) 
#inverse weight by rmse (larger values are smaller weights)
win$rmsewt = 1/win$rmse/sum(1/win$rmse)

win$modnum=1:nrow(win)

#select weight
use.wt='rmsewt'

for(m in seq_along(effects)){
  effects[[m]]$w=win[,use.wt][m]
}

#print weighted mean of windows...
print(
  apply(win[,c('a','p','c')],2,weighted.mean,w=win[,use.wt])
)

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


mean.plt=list()

for(d in names(preds)){

preds[[d]]$m_est =  rowSums(
  as.data.frame(
      lapply(effects,FUN=function(e)
            wtmn(e,d,e$w))
))


  
preds[[d]]$m_down=rowSums(
  as.data.frame(lapply(
  effects,FUN=function(e)
    wtquant(e,d,e$w,0.025))
          )
)

preds[[d]]$m_up=rowSums(
  as.data.frame(lapply(
  effects,FUN=function(e)
    wtquant(e,d,e$w,0.975))
  )
)

#preds$cohort=1:nrow(preds)
#predsm=preds[,! colnames(preds) %in% c('up','down')]
#pp=melt(predsm,id='cohort')

pltrange=range(unlist(
  lapply(preds,function(x) range(x %>% select(-id)))
))

mean.plt[[d]]=ggplot(preds[[d]],
         aes_string(y='m_est',x='id')) + 
         geom_point() + 
         geom_errorbar(aes(ymax=m_up,ymin=m_down)) +
         geom_line(aes(y=actual))

       
}

pdf(file=paste0(imdir,'mean-fit.pdf'))

  grid.arrange(mean.plt[['a']] + ylim(pltrange),
               mean.plt[['p']] + ylim(pltrange),
               mean.plt[['c']] + ylim(pltrange),
               ncol=3)

dev.off()


##posterior for mean effects

#draw list for posterior sample
#post.size=10000
post.size=1000
ytilde = matrix(as.numeric(NA),nrow(tdat),post.size)
post.mods = sample(win$modnum,size=post.size,
                   replace=TRUE,prob=win$rmsewt)

#arror here
for(i in seq_along(post.mods)){
  i.win=win[post.mods[i],]
  modnum=i.win$modnum
  xhat=tdat[,c('a','p','c')]
  #calculate yhat|a single draw from the model
    xhat$a=window(tdat$a,breaks=breaks$a[[modnum]])
    xhat$p=window(tdat$p,breaks=breaks$p[[modnum]])
    xhat$c=window(tdat$c,breaks=breaks$c[[modnum]])

  xhat = as.data.frame(xhat)
  xhat = model.matrix(~.,xhat)
  
  #draw single set of betas and calculate yhat
  #these betas might be good to use for estimating
  #the mean and standard error (instead of the weighted avgs 
  #above)
  b = allmods[[modnum]]$beta[sample(1:1000,1),]
  s2 = allmods[[modnum]]$sigma
  ytilde[,i] = xhat %*% b + rnorm(nrow(xhat),mean=0,sd=s2)
  
  if(i%%100==0){
    cat(i, 'of', length(post.mods),'\n')
  }
  
}

sink(paste0(outdir,'mean-fit-posterior-pval.txt'),type=c('output','message'))

#omnibus bayesian pvalues
print('Summary of acutal to posterior')
print(summary(tdat[,dv]))
print(summary(as.vector(ytilde)))

#mean
print('omnibus bayesian p-value of mean')
sum(apply(ytilde,2,mean)<mean(tdat$y1))/ncol(ytilde)
#sum(apply(ytilde,2,max)<max(tdat[,dv]))/ncol(ytilde)
#sum(apply(ytilde,2,min)>min(tdat$y2))/ncol(ytilde)


#bayesian p-values by period 
library(dplyr)

print('period pvalues')
actual.period = aggregate(tdat[,dv],by=list(tdat$p),mean)
ytilde.period =  sapply(1:ncol(ytilde),function(i)
  aggregate(ytilde[,i],by=list(tdat$p),mean)[[2]]
)
sapply(seq_along(actual.period[[1]]),function(i)
  sum(ytilde.period[i,]>actual.period[i,2])/ncol(ytilde)
)


#bayesian p-values by cohort 
print('cohort pvalues')

actual.cohort = aggregate(tdat[,dv],by=list(tdat$c),mean)
ytilde.cohort =  sapply(1:ncol(ytilde),function(i)
  aggregate(ytilde[,i],by=list(tdat$c),mean)[[2]]
)
sapply(seq_along(actual.cohort[[1]]),function(i)
  sum(ytilde.cohort[i,]>actual.cohort[i,2])/ncol(ytilde)
)


#bayesian p-values by age 
print('age pvalues')

actual.age = aggregate(tdat[,dv],by=list(tdat$a),mean)
ytilde.age =  sapply(1:ncol(ytilde),function(i)
  aggregate(ytilde[,i],by=list(tdat$a),mean)[[2]]
)
sapply(seq_along(actual.age[[1]]),function(i)
  sum(ytilde.age[i,]>actual.age[i,2])/ncol(ytilde)
)
sink()


predy.age = apply(ytilde.age,1,function(x)
  c(mean=mean(x),up=quantile(x,0.975),low=quantile(x,0.025)))
predy.age = as.data.frame(t(predy.age))
predy.age$a = c(1:nrow(predy.age)); colnames(predy.age) = c('mean','up','low','a')

predy.period = apply(ytilde.period,1,function(x)
  c(mean=mean(x),up=quantile(x,0.975),low=quantile(x,0.025)))
predy.period = as.data.frame(t(predy.period))
predy.period$p = c(1:nrow(predy.period)); colnames(predy.period) = c('mean','up','low','p')

predy.cohort = apply(ytilde.cohort,1,function(x)
  c(mean=mean(x),up=quantile(x,0.975),low=quantile(x,0.025)))
predy.cohort = as.data.frame(t(predy.cohort))
predy.cohort$c = c(1:nrow(predy.cohort)); colnames(predy.cohort) = c('mean','up','low','c')


ggplot(predy.age,aes(x=a,y=mean)) + 
  geom_boxplot(data=tdat, aes(x=a,y=y,group=a),alpha=.05) +
  geom_line() + 
  geom_ribbon(aes(ymin=low,ymax=up),alpha=0.4) 

ggplot(predy.period,aes(x=p,y=mean)) + 
  geom_boxplot(data=tdat, aes(x=p,y=y,group=p),alpha=.05) +
  geom_line() + 
  geom_ribbon(aes(ymin=low,ymax=up),alpha=0.4) 
  
ggplot(predy.cohort,aes(x=c,y=mean)) + 
  geom_boxplot(data=tdat, aes(x=c+20,y=y,group=c),alpha=.05) +
  geom_line() + 
  geom_ribbon(aes(ymin=low,ymax=up),alpha=0.4) 


#print(act)

#predict = ggplot(ytilde.age, aes())

#should do this dynamicall...
yl=range(c(actual.age$x,
             actual.period$x,
             actual.cohort$x))

pdf(paste0(imdir,'mean-fit-post.pdf')) #--conditional??
  par(mfrow=c(1,3))
  plot(actual.age,type='l',ylim=yl);
    lines(actual.age[[1]],apply(ytilde.age,1,mean),type='p')
  plot(actual.period,type='l',ylim=yl); 
    lines(actual.period[[1]],apply(ytilde.period,1,mean),type='p')
  plot(actual.cohort,type='l',ylim=yl); 
    lines(actual.cohort[[1]],apply(ytilde.cohort,1,mean),type='p')
dev.off()

####
#testing of final model

fin.bayes=allmods[[length(allmods)]]
fin.freq = lm(tdat[,dv]~model.matrix(~.,x))

print(summary(fin.freq)$r.squared)
print(mean(fin.bayes$r2))

print(logLik(fin.freq))
print(mean(fin.bayes$ll))

print(
  sum(dnorm(y,mean=predict(fin.freq),sd=mean(fin.bayes$sigma),log=TRUE))
)

window.summary =  win %>% select(a,p,c) %>%
    summarize_all(funs(
      m.wn=mean(.),
      rmse.win=weighted.mean(.,w=win$rmsewt),
      bic.win=weighted.mean(.,w=win$wt),
      bic_prime.win=weighted.mean(.,w=win$w_prime))
      ) 
