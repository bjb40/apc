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

#set betas
#t.beta=data.frame(
#  a=0.3,a2=-0.01,p=-0.04,p2=0.02,c=0.35,c2=-0.0015
#)

t.beta=data.frame(
  a=0.5,a2=-0.1,p=0.04,p2=0.02,c=0,c2=0
)

#set r2
r2=0.2

evar=(1-r2)/r2
e=rnorm(n,0,sqrt(evar))

#simulate and save data
dat$y1 = as.matrix(dat)%*%t(t.beta) + e
tdat = dat %>% select(-c,-a2,-p2,-c2)

save(tdat,file=paste0(datdir,'nsim.RData'))
