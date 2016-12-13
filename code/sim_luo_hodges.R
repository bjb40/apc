#Dev R 3.3.0 "Supposedly Educational"
#simulates data for 4 scenarios from Luo & Hodges SMR 2015; p. 14
#
#Bryce Bartlett

#clear cache
rm(list=ls())

#@@
#preliminaries
library(dplyr); library(ggplot2)

#specific source functions (will become package)
source('funs.R')

#@@helper funcitons

#generate combinations of 20 ages and 20 periods
a=rep(1:20,20)
p=a[order(a)]
c=p-a

combos=data.frame(
  a=a,a2=a^2,
  p=p,p2=p^2,c=c,c2=c^2
)

#generate "x" data
xcov=cov(combos[c('a','p','c')])
Omega=cov2cor(xcov)
combos$x.ind=
  round(runif(nrow(combos),max=max(c),min=min(c)))
  combos$x.ind2=combos$x.ind^2
#add a tiny bit of noise that breaks linear dependence
combos$x.cor = c + runif(nrow(combos),max=1,min=-1)
combos$x.cor2=combos$x.cor^2

 #rm(a,p,c)

#repeat data for 25 individuals
dat=combos; for(i in 2:25){dat=rbind(dat,combos)}; #rm(combos,i)
#confirm only 400 unique combos
cat('400 Unique Observations Test:',
  nrow(unique(dat)) == 400, '\n\n'
)

n=nrow(dat)

#id true betas for 4 scenarios (eq. 8-11)

t.beta=data.frame(
  a=0.3,a2=-0.01,p=0.04,p2=0.02,c=0.35,c2=-0.0015
)
t.beta[2,] = t.beta[1,]; t.beta[2,c('a','a2')]=0
t.beta[3,] = t.beta[1,]; t.beta[3,c('p','p2')]=0
t.beta[4,] = t.beta[1,]; t.beta[4,c('c','c2')]=0

print(t.beta)

#draw unique y for each scenario y to y4 in a list of 100 times 
#/innefficeint memory -- should only store ys
dats=list(); inddats=list(); cordats=list()

for(N in 1:100){
  dats[[N]] = dat[,c('a','a2','p','p2','c','c2')]
  inddats[[N]] = dat[,c('a','a2','p','p2','x.ind','x.ind2')]
  cordats[[N]] = dat[,c('a','a2','p','p2','x.cor','x.cor2')]
  for(s in 1:nrow(t.beta)){
    b = as.matrix(t.beta[s,])
    x=as.matrix(dat[,1:length(b)])
    x.ind=as.matrix(inddats[[N]][,1:length(b)])
    x.cor=as.matrix(cordats[[N]][,1:length(b)])
    
   dats[[N]][,paste0('y',s)]=rnorm(n,mean=x %*% t(b),sd=1)
   inddats[[N]][,paste0('y',s)]=rnorm(n,mean=x.ind %*% t(b),sd=1)
   cordats[[N]][,paste0('y',s)]=rnorm(n,mean=x.cor %*% t(b),sd=1)
  }
}

#####
#Save and export simulated data
#####



  
#plot facets for drawn data of four scenarios (loess smoothing)

#a = reshape(dat,idvar=)

d=ggplot(dats[[1]],aes(x=a)) + 
  geom_smooth(aes(y=y1,col='senario 1')) +
  geom_smooth(aes(y=y2,col='senario 2')) +
  geom_smooth(aes(y=y3,col='senario 3')) +
  geom_smooth(aes(y=y4,col='senario 4')) 
  
print(d)

#id time for OLS estimation of single window constraint
t=Sys.time()

#cut window constraints (per luo and hodges)
x = dats[[1]] %>%
  mutate(
    a=factor(a),
    p=cut_interval(p,length=2),
    c=cut_interval(c,length=5)
  ) %>%
  select(a,p,c,y1,y2,y3,y4)

e.beta = list(
  lm(y1~a+p+c,data=x),
  lm(y2~a+p+c,data=x),
  lm(y3~a+p+c,data=x),
  lm(y4~a+p+c,data=x)
  
)

delta.t = Sys.time()-t

cat('One OLS analysis of a single instance window constraint took:', delta.t,'seconds \n\n')
e.time=100*delta.t
cat('Estimated time for one set of constraints:',e.time/60,'minutes \n\n')

#even bins only (with remainder) --- ignoring double cohort//need function... 
configs=10^3 #need to really think about this... windows are adjacent.../if I make unequal
cat('Estimated time for',configs,'window constraints',(e.time*configs)/60/60,'hours\n\n')


#some extra crap re ... compoosition by age (Raudenbush-style)
x$age = as.numeric(levels(x$a)[x$a])

print(
x %>% group_by(p) %>% summarize(mean(age))
)

print(
x %>% group_by(c) %>% summarize(mean(age))
)

x1 = x %>% group_by(c) %>% mutate(., mcage = age - mean(age),cage=mean(age))

x1 %>% group_by(c) %>% summarize(reg=mean(age),mc=mean(mcage))

print(
  summary(lm(y1~a+p+c-1,dat=x1))
)

print(
  summary(lm(y1~p+c*mcage+mcage^2-1,dat=x1))
)

####ind v. corr

est_windows = function(){
  
}

comp.dats = list(uncorrelated=inddats[[1]],correlated=cordats[[1]])
for(d in 1:2){
  colnames(comp.dats[[d]])[5:6]= c('c','c2')
}

c.beta = lapply(comp.dats,FUN=function(x) 
  lm(y1~a+a2+p+p2+c+c2,dat=x,x=TRUE))


