rm(list=ls())

source('config~.R')

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
  dplyr::select(egal,age,year,birthyear,female,race,educ) %>%
  rename(a=age,
         c=birthyear,
         p=year)

dat = dat[complete.cases(dat),]

#load fit object -- should be chains
load(paste0(outdir,'empirical_res2017-8-26.RData'))

use.wt = win$wt

#calculate opval
y.m = mean(dat$egal)
yhat = do.call(cbind,lapply(allmods,function(x) x$yhat))
sigma = do.call(rbind,lapply(allmods,function(x) x$sigma))

samps = sample(1:length(allmods),1000,replace=TRUE,prob=use.wt)
#not quite right, b/c has error in yhat
yt.m = vector(mode='numeric',length=length(samps))
for(s in seq_along(yhat)){
  ss = sample(1:1000,1) #random draw from gibbs
  yt = rnorm(nrow(dat),mean=yhat[,samps[s]],sd=sigma[samps[s],ss]) #ppd draw
  yt.m[s] = mean(yt) #mean
}

opv=sum(yt.m>y.m)/length(samps)

#calculate equaltiy across each dimension for plotting
###HERE!!

calc_equal = function(df){
  
  pv = vector(mode='numeric',length=ncol(df))
  pv[1] = 0 #by construciton per paper
  
  #test equality of effects -- inlcudes BOTH constrained and estimated...
  for(col in 2:ncol(df)){
    d = df[,col] - df[,col-1] #calculate delta
    pv[col] = min(c(sum(d>0),sum(d<0)))/1000*.5

  }
  return(pv)
} # end equal function


eqs = list(a=NULL,p=NULL,c=NULL)

for(n in names(eqs)){
  eqs[[n]] = do.call(rbind,lapply(effects,function(x) calc_equal(x[[n]])))
}  

#probability that diff between this and next == 0 i.e. they are eual
p.eq = lapply(eqs,function(x) apply(x,2,weighted.mean,w=win$wt))

#classic cut off for rejecting null; 5%; i.e. probability that they are equal
strs = p.eq[['p']]>0.05; names(strs) = unique(dat$p[order(dat$p)])

plot(1:length(p.eq[['p']]),p.eq[['p']])
  
#prior calculated "best plot"
b.plt
