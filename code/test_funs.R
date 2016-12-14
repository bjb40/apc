#simulates data for 4 scenarios from Luo & Hodges SMR 2015; p. 14
#
#Bryce Bartlett

#clear cache
rm(list=ls())
source('config~.R')

library(microbenchmark)

#load test data
load(paste0(datdir,'testdat.RData'))

testwindow=FALSE; testinvertwindows=FALSE

if(testwindow){
######
#WINDOW OBJECT TESTS
######

a=tdat$a
c=tdat$p-tdat$a

cw=window(c,winlength=3)

#even windows
w=window(a,winlength=5)
print(
  table(a,w)
)

#uneven windows
w=window(a,winlength=3)
print(
  table(a,w)
)

#breaks
w=window(a,breaks=c(0,5,10,15,20))
print(
  table(a,w)
)

#no breaks or window (should throw error)
window(a)
} #

######
#full factor length reversion tests
#this is key in estimating "single" window estimates
######
if(testinvertwindows){

a = tdat$a
w = window(a,winlength=5)

#should throw error:
scopedummy(a)

invert.window=scopedummy(w)

print(
  table(invert.window)
)
}

######
#ols linear model tests
######

tst=apc_lm(y2~a+p+cohort,data=tdat,a='a',p='p')

##need to incorprate these into the tests below

b=coefficients(tst$results)
  c.b = b [grepl('cohort|Intercept',names(b))]
cov=vcov(tst$results)
  c.cov = cov[grepl('cohort|Intercept',rownames(cov)),
      grepl('cohort|Intercept',colnames(cov))]

  

predat=model.matrix(~.,data=as.data.frame(tst$blockdat$c))

predict=predat %*% c.b

#this is _assuming_ at default (left out) variables
#should re-estimate based on means??

v = predat %*% c.cov %*% t(predat)

preds = data.frame(est=predict,se=sqrt(diag(v)))
preds$up = preds$est+preds$se*1.96
preds$down = preds$est-preds$se*1.96

####
#plots -- not necessarily for functions
####

#cohort plots

library(ggplot2)
g=ggplot(preds,
       aes(y=est,x=1:nrow(preds))) + 
       geom_point() + 
       geom_errorbar(aes(ymax=up,ymin=down))

print(g)
