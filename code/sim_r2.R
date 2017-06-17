#Bryce Bartlett
#simulates data for use in apc models

#clear cache
rm(list=ls())
source('config~.R')

#preliminaries
library(dplyr)

a=rep(1:20,20)
p=a[order(a)]
c=p-a

combos=data.frame(
  a=a,a2=a^2,
  p=p,p2=p^2,c=c,c2=c^2
)

dat=combos; for(i in 2:25){dat=rbind(dat,combos)}; #rm(combos,i)
#confirm only 400 unique combos
cat('400 Unique Observations Test:',
    nrow(unique(dat)) == 400, '\n\n'
)

n=nrow(dat)

#set betas (luo)
#t.beta=data.frame(
#  a=0.3,a2=-0.01,p=-0.04,p2=0.02,c=0.35,c2=-0.0015
#)

#t.beta=data.frame(
#  a=0.3,a2=-0.01,p=-0.04,p2=0.02,c=0.35,c2=-0.0015
#)

t.beta=data.frame(t(runif(6,-1,1)))
colnames(t.beta)=c('a','a2','p','p2','c','c2')
#rescale quadratic effects to be smaller
t.beta[,c('a2','p2','c2')] = t.beta[,c('a2','p2','c2')]/10

save(t.beta,file=paste0(datdir,'sim2_tbeta.RData'))

#kill effect
#t.beta[,c('c','c2')] = 0

#set r2--need to fix

r2=0.03

evar=(1-r2)/r2
e=rnorm(n,0,sqrt(evar))

#simulate and save data
dat$s1 = as.matrix(dat)%*%t(t.beta)[,1]
dat$y1 = dat$s1 + e
#tdat = dat
tdat = dat %>% select(-c,-a2,-p2,-c2)

print(
  summary(lm(y1~a+p,data=tdat))
)

save(tdat,file=paste0(datdir,'nsim.RData'))

pltdat=list()

#pull unique data -- should make this a loop
uns = list(
  c = unique(dat$c), a=unique(dat$a), p=unique(dat$p)
)

#need to hold means for margins (i.e. margins command in stata)??
dims=c('a','p','c')

for(d in dims){
  o.dims = dims[which(!dims==d)]
  o.means=apply(dat[,o.dims],2,mean)
  xdat=data.frame(uns[[d]],o.means[1],o.means[2]); colnames(xdat)=c(d,o.dims)
  xdat$a2 = xdat$a^2; xdat$p2 = xdat$p^2; xdat$c2 = xdat$c^2
  x = as.matrix(xdat)
  b = as.matrix(t.beta)[,colnames(x)]
  
  #test for matching
  if(!all(colnames(x) == colnames(b))){stop('Problem with matrix ordering.')}
  
  pltdat[[d]] = data.frame(
    id = uns[[d]],
    s1 = (x %*% b)[,1] 
  )
  
}

save(pltdat,file=paste0(datdir,'nsim_fits.RData'))

#run algorithm
source('algorithm.R')

