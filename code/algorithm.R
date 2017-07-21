###

source('config~.R')

#can pass these to ... parallel chains, maybe
load(paste0(datdir,'nsim.RData'))
tdat$c=tdat$p-tdat$a
dv='y1'
actual='s1'
#y=tdat[,dv]

#begin chain fuction
draw_chains = function(...){

source('config~.R')

#prelim
dv='y1'
actual='s1'

#load test data
#load(paste0(datdir,'testdat.RData')) #this is luo
load(paste0(datdir,'nsim.RData'))
tdat$c=tdat$p-tdat$a

#@@@@@@@@@@@
#sampling from model subspace
#@@@@@@@@@@@

y=tdat[,dv]
maxc= max(tdat$c)
minc = min(tdat$c)

#perturb
"
tdat$c = tdat$c + round(runif(nrow(tdat),-1,1))
tdat$c = ifelse(tdat$c>maxc,maxc,tdat$c)
tdat$c = ifelse(tdat$c<minc,minc,tdat$c)
print('Check whether perturbation solves exact multicolinearity')
print(summary(lm(y1~a+p+c,data=tdat)))
Sys.sleep(5)
"

allmods=list() #may run into size constraints/may need to limit to best mods... 
effects=xhats=ppd=list()
tm=Sys.time()
avtm=0

window.sample=function(var,alph){
  #input is continuous of a,p,c
  #alpha is a vector the length of unique entries in var that is fed to dirichelet
  #output is factor with uniform, random window constraints
  #see dirichelet, e.g. for alternative algorithms
  vals=unique(var) #center on 0; assume continuous
  len=length(vals)
  
  alph=unlist(alph)

  dp=rdirichlet(1,alph)
  #segment the stick
  segments=round(cumsum(dp*len))
  #identify changes in segments
  sb=segments
  for(i in seq_along(segments)){
    if(i==1){next}
    sb[i] = segments[i]-segments[i-1]
  }
  
  #id inclusive breaks
  breaks=vals[sb>=1]
  
  #because breaks are *inclusive*, must include min-1; ensure max
  if(min(breaks)>(min(var)-1)){
    breaks=c((min(var)-1),breaks)
  }
  if(max(breaks)<max(var)){
    breaks=c(breaks,max(var))
  }
  
  wins=window(var,breaks=breaks)

  return(wins)
}

#set of numbers of random samples
n.samples=100

##you are using the wrong test --- for MC3, should be bic approx to bayes factor
#see raferty
#also, your proposal for alpha can use a mulitvariate...


#holder df for model summary data 
win = data.frame(a=numeric(), p=numeric(), c=numeric())
breaks=list(a=list(),p=list(),c=list())

d = c('a','p','c')
dl = c(length(unique(tdat$a)),length(unique(tdat$p)),length(unique(tdat$c)))
names(dl) = d

#set starting values
all.alphas = lapply(d,function(x) data.frame(t(rep(dl[x],length(unique(tdat[,x]))))))
#all.nwins = lapply(d,function(x) 3)
#all.nwins = lapply(d,function(x) length(unique(tdat[,x]))/2)
#all.nwins = lapply(d,function(x) length(unique(tdat[,x]))-1)

names(all.alphas) = d #names(all.nwins) = d

#accept rate
acc=0
#count boundary conditions rejections
bound=0

#mc3 sampler (prior model probabilities are equal)
for(s in 2:n.samples){

  #reset dataframe
  x=tdat[,c('a','p','c')]
  
  #draw from proposal distributions
  #all.nwins = lapply(all.nwins, function(x)
  #                        append(x,x[s-1]+rnorm(1,mean=0,sd=1)))
  #all.nwins = lapply(list(a='a',p='p',c='c'),function(x)
  #  append(all.nwins[[x]],runif(1,2,length(unique(tdat[,paste(x)]))))
  #  )

  
  all.alphas= lapply(all.alphas, function(x)
                    rbind(x,x[s-1,]+rnorm(ncol(x),mean=0,sd=1)))
  
  #all.alphas= lapply(all.alphas, function(x)
  #  rbind(x,runif(ncol(x),0,10)))

  for(d in seq_along(all.alphas)){rownames(all.alphas[[d]]) = 1:nrow(all.alphas[[d]])}  

  #if(any(unlist(all.alphas)<0 | any(unlist(all.nwins)<2))){
  if(any(unlist(all.alphas)<0)){
      
    bound=bound+1
    out.al=sum(unlist(all.alphas)<0)
    #out.wi=sum(unlist(all.nwins)<2)
    #print(c(out.al,out.wi))
    #s=s-1
    #mnum=mnum-1 #should consolidate these
    cat('\n\nOut-of-Sample-Space Warning.\n\n')
    #acc=acc-1
    for(d in seq_along(all.alphas)){
      #all.nwins[[d]][s]=all.nwins[[d]][s-1]
      all.alphas[[d]][s,]=all.alphas[[d]][s-1,]
     #note that this samples different windows with same hyper param
    }

  }
  
  #draw random window samples
  x$a=window.sample(x$a,all.alphas$a[s,]) 
  x$p=window.sample(x$p,all.alphas$p[s,]) 
  x$c=window.sample(x$c,all.alphas$c[s,]) 
  
  #skip if unideintified
  la = length(levels(x$a)) == length(unique(tdat$a))
  lp = length(levels(x$p)) == length(unique(tdat$p))
  lc = length(levels(x$c)) == length(unique(tdat$c))
  if(all(la,lp,lc)){
    for(d in seq_along(all.alphas)){
      #all.nwins[[d]][s]=all.nwins[[d]][s-1]
      all.alphas[[d]][s,]=all.alphas[[d]][s-1,]
    }
  }
  
  #collect model data
  nr=data.frame(a=length(levels(x$a)),
                p=length(levels(x$p)),
                c=length(levels(x$c)))
  win=rbind(win,nr)
  
  #collect breaks data
  breaks$a[[s]]=attr(x$a,'breaks')
  breaks$p[[s]]=attr(x$p,'breaks')
  breaks$c[[s]]=attr(x$c,'breaks')

"  
if(s%%10==0){
      cat('\n\nEstimating model',length(allmods),
            'with windows:',
            '\n\tage    ',nr$a,
            '\n\tperiod ',nr$p,
            '\n\tcohort ',nr$c,
            '\n',
          '\nwindows:\n\t',
          unlist(lapply(all.nwins, function(x) round(x[s],2))),
          '\nalphas\n\tmin:\n',
              '\t',
  unlist(lapply(all.alphas,function(x) round(min(unlist(x[s,])),4))),
  '\n\tmax:\n',
  '\t',
  unlist(lapply(all.alphas,function(x) round(max(unlist(x[s,])),4))),
  
        '\nacceptance rate:',(acc/s),
        '\nboundary rate:', (bound/s),
        '\nAverage model time:',avtm,
        '\n\n')
}
"
      #add esitmate / lin_gibbs is bayesian gibbs sampler
      mnum = length(allmods)+1
      
      #reassign random references to each vector
      a.lev=length(levels(x$a)); a.b = sample(1:a.lev,1)
      p.lev=length(levels(x$p)); p.b = sample(1:p.lev,1)
      c.lev=length(levels(x$c)); c.b = sample(1:c.lev,1)
      
      form.c = as.formula(paste("~C(a,contr.treatment(a.lev,base=a.b))+
        C(p,contr.treatment(p.lev,base=p.b))+
        C(c,contr.treatment(c.lev,base=c.b))"))
      
      xmat = model.matrix(form.c,data=x)
      m = allmods[[s]] = lin_gibbs(y=y,x=xmat)
      
      #this seems to match the 's1' margin property
      #i don't acutally have to do the scopedummy stuff, b/c I can use the yhat
      #from the lin_gibbs/same with the sampling below
      grand.means = data.frame(
        a = window(mean(tdat$a),breaks=attr(x$a,'breaks')),
        p = window(mean(tdat$p),breaks=attr(x$p,'breaks')),
        c = window(mean(tdat$c),breaks=attr(x$c,'breaks'))
      )
      
      grand.means$a=relevel(grand.means$a,ref=a.b)
      grand.means$p=relevel(grand.means$p,ref=p.b)
      grand.means$c=relevel(grand.means$c,ref=c.b)
      
      grand.means=(model.matrix(form.c,grand.means))
      
     
      blockdat=lapply(x,scopedummy)
      blockdat$a = relevel(blockdat$a,ref=a.b)
      blockdat$p = relevel(blockdat$p,ref=p.b)
      blockdat$c = relevel(blockdat$c,ref=c.b)
      
      predat=lapply(blockdat,FUN=function(x) 
        model.matrix(~.,data=as.data.frame(x)))

      effects[[mnum]] = xhats[[mnum]] = list()
      #xhat = list()

      betas=list()
      for(eff in names(predat)){
        #fix colnames
        colnames(predat[[eff]]) = sub('x',eff,colnames(predat[[eff]]))
        
        #calculate means for xhat, & id effects at issue---this was replaced...
        xhat=grand.means[rep(seq(nrow(grand.means)), nrow(predat[[eff]])),]

        #replace means of effect dimensions with indicator in matrix
        calceff = grepl(paste0(eff,'.lev|Intercept'),colnames(xhat))
        xhat[,calceff] = predat[[eff]]
        xhats[[mnum]][[eff]] = xhat
  
        
        effects[[mnum]][[eff]] = t(xhat %*% t(m$betas))
        colnames(effects[[mnum]][[eff]]) = paste0(eff,unique(tdat[,eff]))
      }
      
      effects[[mnum]]$bic=m$bic
      
        avtm=(avtm*(length(allmods)-1)+Sys.time()-tm)/length(allmods)
        tm=Sys.time()
      
        
        if(s==2){next}
        
        #selection criterion
        
        #bayes factorf approximation
        bf=exp((allmods[[s]]$bic-allmods[[s-1]]$bic)/2)
        R = min(1,bf)
        #R = min((allmods[[s]]$bic/ allmods[[s-1]]$bic),1)
        #print(R)
        if (R < runif(1)){
          acc = acc+1
        } else {
          
          for(d in seq_along(all.alphas)){
            #all.nwins[[d]][s]=all.nwins[[d]][s-1]
            all.alphas[[d]][s,]=all.alphas[[d]][s-1,]
          }

          
        }
        
        
}#end sampling loop


allmods = allmods[2:length(allmods)]
effects = effects[2:length(effects)]
xhats = xhats[2:length(xhats)]
breaks = lapply(breaks, function(x)
  x[2:length(x)])


res = list(
  allmods=allmods,
  effects=effects,
  xhats=xhats,
  breaks=breaks,
  win=win,
  n.samples=125,
  acc=acc,
  bound=bound
)

return(list(res))

}

library(parallel)

#draw 4 chains

print('Begin drawing chains...')

cl <- makeCluster(mc <- getOption("cl.cores", 4))
chains=do.call(base::c,parLapply(cl=cl,1:4,draw_chains))

#end cluster
stopCluster(cl)

#combine chains --- can calculate R-hats here...
extract = function(l,name,as.df=FALSE){
  #extracts and combines name object from list 'l'
  res=do.call(base::c,lapply(l,function(x) x[[name]]))
  if(as.df){
    res=do.call(rbind,lapply(l,function(x) x[[name]]))
  }
  
  return(res)
}

allmods = extract(chains,'allmods')
effects = extract(chains,'effects')
xhats = extract(chains,'xhats')
breaks = extract(chains,'breaks')
win = extract(chains,'win',as.df=TRUE)
n.samples = sum(unlist(extract(chains,'n.samples')))
acc = sum(unlist(extract(chains,'acc')))
bound = sum(unlist(extract(chains,'bound')))
rm(chains)

##post-processing --- best model
bics=unlist(lapply(allmods,FUN=function(x) x$bic))
bics_prime=unlist(lapply(allmods,FUN=function(x) x$bic_prime))
r2=unlist(lapply(allmods,FUN=function(x) mean(x$r2)))

best=which(bics_prime==min(bics_prime))

library(ggplot2)
library(gridExtra)

#load(paste0(datdir,'luo_sim_fits.RData'))
load(paste0(datdir,'nsim_fits.RData'))

best.plt = list()
preds = list()

#for(d in c('a','p','c')){
for(d in c('a','p','c')){
  if(length(effects[[best]][[d]])==0){
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
  preds[[d]]$actual=pltdat[[d]]$s1[order(pltdat[[d]]$id)]
  
  preds[[d]]$id=unique(tdat[,d])[order(pltdat[[d]]$id)]
  #preds[[d]]=preds[[d]] %>% arrange(id)
  
  
}


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
win$mse = unlist(lapply(allmods,FUN=function(x) mean(x$rmse^2))) 
win$msewt=1/win$mse/sum(1/win$mse)

win$modnum=1:nrow(win)

#select weight
#use.wt='rmsewt'
use.wt='wt'

#generate sims data if it doesn't exist
#place summary, including pvalues, summary stats, estimates, and true beta...
srf = paste0(outdir,'simres.RData')
if(!file.exists(srf)){
  simres = list()
  save(simres,file=srf)
  simnum = 1
} else{
  load(srf)
  simnum = length(simres)+1
}

#incomplete list... others self-explanatory
tempsim = list(tbeta = NULL, #true beta
               r2 = NULL, #summary of r2 values
               bic = NULL, #summar of BIC values
               pvals = NULL, #list of bayesian p-values --- omnibus + sliced
               wt = NULL, #weighting scheme
               preds = NULL, #mean and best estimates w/actual and overlap (actual can be adj)
               ppdsum = NULL, #summary of ppd 
               datsum = NULL #summary of data
)

tempsim$wt = use.wt
tempsim$samp = n.samples
tempsim$acc = acc/n.samples
tempsim$bound = bound/n.samples

for(mod in seq_along(effects)){
  effects[[mod]]$w=win[,use.wt][mod]
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


       
}

#now = Sys.time()

cl <- makeCluster(mc <- getOption("cl.cores", 4))
clusterExport(cl=cl, varlist=c("win", "tdat","allmods",'use.wt'))

draw_post=function(...){
#draw list for posterior sample
  post.size=250
  
  ytilde = matrix(as.numeric(NA),nrow(tdat),post.size)
  post.mods = base::sample(win$modnum,size=post.size,
                   replace=TRUE,prob=win[,use.wt])
  
  #drawing posterior samples (with replacement)
  for(i in seq_along(post.mods)){
    i.win=win[post.mods[i],]
    modnum=i.win$modnum
  
    ytilde[,i] = allmods[[modnum]]$yhat + 
      rnorm(nrow(ytilde),mean=0,sd=allmods[[modnum]]$sigma)
  }
  
  return(list(ytilde))
}

print('Begin drawing posterior samples...')

post_chains=do.call(base::c,parLapply(cl=cl,1:4,draw_post))

#end cluster
stopCluster(cl)

#print(Sys.time()-now)

ytilde = do.call(cbind,post_chains)
rm(post_chains)

print('Finalize post-processing, and output results...')

sink(paste0(outdir,'mean-fit-posterior-pval.txt'),type=c('output','message'))

#omnibus bayesian pvalues
#print('Summary of acutal to posterior')
tempsim$datsum = summary(tdat[,dv])
tempsim$ppdsum = summary(as.vector(ytilde))
print(tempsim$datsum)
print(tempsim$ppdsum)

tempsim$pvals = list()

#mean
print('omnibus bayesian p-value of mean')
tempsim$pvals[['omnibus']] = sum(apply(ytilde,2,mean)<mean(tdat$y1))/ncol(ytilde)
print(tempsim$pvals[['omnibus']])

sink()

#sum(apply(ytilde,2,max)<max(tdat[,dv]))/ncol(ytilde)
#sum(apply(ytilde,2,min)>min(tdat$y2))/ncol(ytilde)

#bayesian p-values by period 

#print('period pvalues')
actual.period = aggregate(tdat[,dv],by=list(tdat$p),mean)
ytilde.period =  sapply(1:ncol(ytilde),function(i)
  aggregate(ytilde[,i],by=list(tdat$p),mean)[[2]]
)

tempsim$pvals[['p']] = 
  sapply(seq_along(actual.period[[1]]),function(i)
    sum(ytilde.period[i,]>actual.period[i,2])/ncol(ytilde)
  )


#bayesian p-values by cohort 
#print('cohort pvalues')

actual.cohort = aggregate(tdat[,dv],by=list(tdat$c),mean)
ytilde.cohort =  sapply(1:ncol(ytilde),function(i)
  aggregate(ytilde[,i],by=list(tdat$c),mean)[[2]]
)
tempsim$pvals[['c']] =
  sapply(seq_along(actual.cohort[[1]]),function(i)
    sum(ytilde.cohort[i,]>actual.cohort[i,2])/ncol(ytilde)
  )


#bayesian p-values by age 
#print('age pvalues')
actual.age = aggregate(tdat[,dv],by=list(tdat$a),mean)
ytilde.age =  sapply(1:ncol(ytilde),function(i)
  aggregate(ytilde[,i],by=list(tdat$a),mean)[[2]]
)
tempsim$pvals[['a']] =
  sapply(seq_along(actual.age[[1]]),function(i)
    sum(ytilde.age[i,]>actual.age[i,2])/ncol(ytilde)
  )

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


preda = ggplot(predy.age,aes(x=a,y=mean)) + 
  geom_boxplot(data=tdat, aes(x=a,y=tdat[,dv],group=a),alpha=.05) +
  geom_line() + 
  geom_ribbon(aes(ymin=low,ymax=up),alpha=0.4) 

#Sys.sleep(5)

predp = ggplot(predy.period,aes(x=p,y=mean)) + 
  geom_boxplot(data=tdat, aes(x=p,y=tdat[,dv],group=p),alpha=.05) +
  geom_line() + 
  geom_ribbon(aes(ymin=low,ymax=up),alpha=0.4) 

#Sys.sleep(5)
  
predc = ggplot(predy.cohort,aes(x=c,y=mean)) + 
  geom_boxplot(data=tdat, aes(x=c+20,y=tdat[,dv],group=c),alpha=.05) +
  geom_line() + 
  geom_ribbon(aes(ymin=low,ymax=up),alpha=0.4) 

pdf(file=paste0(imdir,'comp-post.pdf'))

grid.arrange(preda + theme_minimal(),
             predp + theme_minimal(),
             predc + theme_minimal(),
             nrow=3)

dev.off()

tempsim$tbeta = t.beta
tempsim$preds = preds
tempsim$r2 = summary(r2)
tempsim$bic = summary(bics)

#save summary of simulation data
simres[[simnum]] = tempsim
save(simres,file=srf)