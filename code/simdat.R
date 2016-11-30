#Dev R 3.3.0 "Supposedly Educational"
#simulates age-period-cohort effects on continuous outcome; ols
#
#Bryce Bartlett

#clear cache
rm(list=ls())

#@@@@@@@@
#funcitons and preliminaries
outdir="H:/projects/APC/output/"

#@@@@@@@@
#dependencies
library(lme4); library(ggplot2); library(dplyr)


#@@@@@@@@
#id variance and effects

#need to use Luo and Hodges 2015 estimates
beta=c(.5,-.7,0,.6); names(beta)=c('int','age','period','cohort')
#set r-squared for simulation
#r2=.1
sig=.15
n=10000 #sample size

#sumulate data
simdat = data.frame(
  'int'=rep(1,n),
  'age'=rpois(n,20),
  'period'= sample(0:30,n,replace=TRUE)
)

simdat$cohort = simdat$period-simdat$age

#draw outcome
simdat$y=rnorm(
  n=nrow(simdat),mean=as.matrix(simdat) %*% beta,
sd=sig)

print(summary(simdat))

print(hist(simdat$y))

print(
summary(lm(y~age+period+cohort,data=simdat))
)

print(beta)

#mean center period and cohrt
simdat$c.cohort = mean(simdat$cohort) -  simdat$cohort
simdat$c.period = mean(simdat$period) -  simdat$period

#HAPC
hapc = lmer(y~-1+age + c.cohort + c.period + (1 | c.cohort) + (1 | c.period),simdat)
summary(hapc)

#@@@@@@@@@@@@@@@
#windows function
#@@@@@@@@@@@@@@@

#helper function for slicing
window=function(var,w){
  #input is a numeric vector for a,p,c
  #w is a positive integer for the size of the "window" for dummy variables
  #output is an ordered factor with the windows...
  
  #calculate nubmer of "breaks" required
  u=length(unique(var)) #uniqu items
  b=round(u/w)
  
  return(cut(var,b))
  
}

ndat = data.frame(y=simdat$y,
                  age=simdat$age,
                  period=window(simdat$period,3),
                  cohort=window(simdat$cohort,3))

ndat$period=as.numeric(ndat$period); ndat$cohort = as.numeric(ndat$cohort)

#age2 for linear example
simdat$age2=simdat$age-18
  
print(
  summary(lm(y~age+period+cohort,data=ndat))
)

write.csv(simdat,file=paste0(outdir,'simdat.csv'))


#restrictions test
#@@@@

#if only linear, then, where p=a+c
#yhat1=ba*a + bc*c 
#yhat2=bp*p = bp*(a+c)
#yhat1-yhat2=ba*a + bc*c - bp*a - bp*c = 0  ??

simdat$mcy = simdat$y-mean(simdat$y) # doesn't work

simdat$invage=simdat$age*-1
m1=lm(y~age+cohort,data=simdat)
m2=lm(y~period,data=simdat)
m3=lm(y~period+invage,data=simdat)

b1=coef(m1)
v1=vcov(m1)
b2=coef(m2)
b2a=c(b2,b2[length(b2)])
v2=vcov(m2)

b3=coef(m3)
v3=vcov(m3)

x1=model.matrix(~age+cohort,data=simdat)
x2=model.matrix(~period,data=simdat)

yhat1=x1%*%b1
#vhat1=x1%*%v1%*%t(x1)

yhat2=x2%*%b2
yhat2a=x1%*%b2a

print(
  mean(yhat1-yhat2)
)


p3=predict.lm(m3,se.fit=TRUE)

#yhat=xhat%*%b

#v=xhat%*%v%*%t(xhat)

