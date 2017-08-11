#Load previously subset gss for gender egal
#
#Bryce Bartlett

rm(list=ls())

source('config~.R')

#previously cleaned data for gender egalitarian project
dat.f = "H:/projects/gender_egal/output/private~/subdat.csv"

if(file.exists(dat.f)){
  dat = read.csv(dat.f)
  #save object for runing on the server
  save(dat,file=paste0(outdir,'~dat.RData')) #~ keeps it from pushing to git
} else{
  load(paste0(outdir,'~dat.RData'))
}

dat = dat %>%
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
  dplyr::select(egal,age,year,birthyear,female,race,educ)

t=nrow(dat)

dat = dat[complete.cases(dat),]

cat('Proportion non-missing:',nrow(dat)/t)


library(parallel)


draw_chains = function(...){
   source('config~.R')

  y = dat$egal
  #y = dat$happy
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
n.samples=250

#holder df for model summary data 
win = data.frame(a=numeric(), p=numeric(), c=numeric())
breaks=list(a=list(),p=list(),c=list())

d = c('a','p','c')
dl = c(length(unique(dat$a)),length(unique(dat$p)),length(unique(dat$c)))
names(dl) = d

#set starting values
all.alphas = lapply(d,function(x) data.frame(t(rep(dl[x]/2,length(unique(dat[,x]))))))

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
    rbind(x,x[s-1,]+rnorm(ncol(x),mean=0,sd=0.5)))

  for(d in seq_along(all.alphas)){rownames(all.alphas[[d]]) = 1:nrow(all.alphas[[d]])}  
  
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
  #blockdat$p = scopedummy(w=x$p,unique.vals=unique(dat$p))
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
  
  #bayes factor approximation
  bf=exp((allmods[[s]]$bic-allmods[[s-1]]$bic)/2)
  R = min(1,bf)
  #print(R)
  if (R < runif(1)){
    acc = acc+1
  } else {
    
    for(d in seq_along(all.alphas)){
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

return(res)

} #end draw chains



#define cluster to draw samples
cat('\n\nBegin drawing chains...')

cl <- makeCluster(mc <- getOption("cl.cores", 4))
clusterExport(cl=cl, varlist=c("dat"))
chains=do.call(list,parLapply(cl=cl,1:4,draw_chains))

#end cluster
cat('\n\nEnd drawing chains...')
stopCluster(cl)

save.image(file=paste0(outdir,'empirical_res.RData'))


