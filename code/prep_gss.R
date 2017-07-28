#Load previously subset gss for gender egal
#
#Bryce Bartlett

rm(list=ls())

source('config~.R')

#previously cleaned data for gender egalitarian project
dat = read.csv("H:/projects/gender_egal/output/private~/subdat.csv") %>%
  mutate(happy=3-happy,
         fechld = 4-fechld,
         fefam = fefam-1,
         fepol = fepol-1,
         fepresch = fepresch-1,
         egal = fechld + fepresch + fefam + fepol,
         birthyear = year - age,
         female=sex-1,
         race=factor(race,labels=c('White','Black','Other')))

#limit to complete cases // inlcude educaton only as covariate
dat = dat %>%
  dplyr::select(happy,age,year,birthyear,female,race,educ)

t=nrow(dat)

dat = dat[complete.cases(dat),]

remain = t- nrow(dat)
cat('Proportion non-missing:',remain/t)

#y = dat$egal
y = dat$happy
allmods=list() #may run into size constraints/may need to limit to best mods... 
effects=xhats=ppd=list()
tm=Sys.time()
avtm=0

#draw_chains = function(...){
#   source('config~.R')
  

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
  #print(breaks)
  
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

dat = dat %>%
  rename(a=age,
         p=year,
         c=birthyear)

#set of numbers of random samples
n.samples=150

#holder df for model summary data 
win = data.frame(a=numeric(), p=numeric(), c=numeric())
breaks=list(a=list(),p=list(),c=list())

d = c('a','p','c')
dl = c(length(unique(dat$a)),length(unique(dat$p)),length(unique(dat$c)))
names(dl) = d

#set starting values
all.alphas = lapply(d,function(x) data.frame(t(rep(dl[x]/2,length(unique(dat[,x]))))))
#all.alphas = lapply(d,function(x) data.frame(rep(1,length(unique(dat[,x])))))


names(all.alphas) = d #names(all.nwins) = d

#accept rate
acc=0
#count boundary conditions rejections
bound=0

#mcmc sampler (prior model probabilities are equal)
for(s in 2:n.samples){
  
  #reset dataframe
  x=dat[,c('a','p','c')]

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
  
  #skip if unideintified
  la = length(levels(x$a)) == length(unique(dat$a))
  lp = length(levels(x$p)) == length(unique(dat$p))
  lc = length(levels(x$c)) == length(unique(dat$c))
  if(all(la,lp,lc)){
    for(d in seq_along(all.alphas)){
      #all.nwins[[d]][s]=all.nwins[[d]][s-1]
      all.alphas[[d]][s,]=all.alphas[[d]][s-1,]
    }
  }
  
  
  #draw random window samples
  x$a=window.sample(x$a,all.alphas$a[s,]) 
  x$p=window.sample(x$p,all.alphas$p[s,]) 
  x$c=window.sample(x$c,all.alphas$c[s,]) 
  
  
  #collect model data
  nr=data.frame(a=length(levels(x$a)),
                p=length(levels(x$p)),
                c=length(levels(x$c)))
  win=rbind(win,nr)
  
  #collect breaks data
  breaks$a[[s]]=attr(x$a,'breaks')
  breaks$p[[s]]=attr(x$p,'breaks')
  breaks$c[[s]]=attr(x$c,'breaks')
  
  #add esitmate / lin_gibbs is bayesian gibbs sampler
  mnum = s
  
  #reassign random references to each vector
  a.lev=length(levels(x$a)); a.b = sample(1:a.lev,1)
  p.lev=length(levels(x$p)); p.b = sample(1:p.lev,1)
  c.lev=length(levels(x$c)); c.b = sample(1:c.lev,1)
  
  form.c = as.formula(paste("~C(a,contr.treatment(a.lev,base=a.b))+
                            C(p,contr.treatment(p.lev,base=p.b))+
                            C(c,contr.treatment(c.lev,base=c.b))"))
  
  xmat = model.matrix(form.c,data=x)
  
  m = allmods[[s]] = lin_gibbs(y=y,x=xmat)
  #m = allmods[[s]] = tryCatch({lin_gibbs(y=y,x=xmat)},
  #                            finally=next)
  
  grand.means = data.frame(
    a = window(mean(dat$a),breaks=attr(x$a,'breaks')),
    p = window(mean(dat$p),breaks=attr(x$p,'breaks')),
    c = window(mean(dat$c),breaks=attr(x$c,'breaks'))
  )
  
  grand.means$a=relevel(grand.means$a,ref=a.b)
  grand.means$p=relevel(grand.means$p,ref=p.b)
  grand.means$c=relevel(grand.means$c,ref=c.b)
  
  grand.means=(model.matrix(form.c,grand.means))
  
  
  blockdat=lapply(x,scopedummy)
  blockdat$p = scopedummy(w=x$p,unique.vals=unique(dat$p))
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
    #Error here??
    #colnames(effects[[mnum]][[eff]]) = paste0(eff,unique(dat[,eff]))
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


#allmods = allmods[2:length(allmods)]
#effects = effects[2:length(effects)]
#xhats = xhats[2:length(xhats)]
#breaks = lapply(breaks, function(x)
#  x[2:length(x)])
breaks = lapply(breaks, function(x)
    x[1:length(x)])

res = list(
  allmods=allmods,
  effects=effects,
  xhats=xhats,
  breaks=breaks,
  win=win,
  n.samples=n.samples,
  acc=acc,
  bound=bound
)

#return(res)

#} #end draw chains

allmods = allmods[2:length(allmods)]

bics=unlist(lapply(allmods,FUN=function(x) x$bic))
bics_prime=unlist(lapply(allmods,FUN=function(x) x$bic_prime))
r2=unlist(lapply(allmods,FUN=function(x) mean(x$r2)))

#NEED TO FIND WHY ERROR AND FIX
best=which(bics_prime==min(bics_prime))
if(is.null(effects[[best]])){best=3}


library(ggplot2)
library(gridExtra)

best.plt = list()
preds = list()

for(d in c('a','p','c')){
#for(d in c('a','p','c')){
#  if(length(effects[[best]][[d]])==0){
#    effects[[best]][[d]]=matrix(0,1000,20)
#    if(d == 'c'){
#      effects[[best]][[d]]=cbind(effects[[best]][[d]],matrix(0,1000,19))
#    }
#  }
  
  preds[[d]]=as.data.frame(apply(effects[[best]][[d]],2,mean))
  colnames(preds[[d]])='est'
  rng=apply(effects[[best]][[d]],2,quantile,c(0.025,0.975))
  preds[[d]]$up = rng[2,]
  preds[[d]]$down = rng[1,]
  #s1-s4 are for scenarios --- needs to match with y1-y4

  #pldat doesn't exist -- need to figure about ordering
  preds[[d]]$id=unique(dat[,d])[order(unique(dat[,d]))]
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

effects = effects[2:length(effects)]

#select weight
#use.wt='rmsewt'
use.wt='wt'

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


df = do.call(rbind,preds)
df$type=substr(rownames(df),1,1)

b.plt=ggplot(df,aes(x=id,y=est)) +
  geom_line(alpha=0.5) +
  geom_ribbon(aes(ymin=down,ymax=up), alpha=0.5)  +
  facet_grid(.~type,scales='free_x') +
  theme_minimal() +
  ylab('Effect') +
  xlab('APC Value')

print(b.plt)

m.plt=ggplot(df,aes(x=id,y=m_est)) +
  geom_line() +
  geom_ribbon(aes(ymin=m_down,ymax=m_up), alpha=0.25) +
  facet_grid(.~type,scales='free_x') +
  theme_minimal() +
  ylab('Effect') +
  xlab('APC Value')

print(m.plt)

save(allmods,
     all.alphas,
     effects,win,
     b.plt,
     n.samples,
     bound,
     acc,
     file=paste0(outdir,'empirical_res.RData'))

