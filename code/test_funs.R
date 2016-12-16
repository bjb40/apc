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

winlist = list(age=4,period=5,cohort=4)

#prepare continuous cohort variable 
#for ineraction (piecewise slopes) -- doesn't work
#tdat$c=tdat$p-tdat$a

tst=apc_lm(y2~a+p+cohort,data=tdat,a='a',p='p')



####
#plots -- not necessarily for functions
####

#cohort plots

library(ggplot2)

load(paste0(datdir,'luo_sim_fits.RData'))

  preds=tst$block.eff[['cohort']]
  preds$up=preds$est + 1.96*preds$se
  preds$down=preds$est - 1.96*preds$se
  preds$actual=pltdat$s1[order(pltdat$cohort)]
  
g=ggplot(preds,
       aes(y=est,x=1:nrow(preds))) + 
       geom_point() + 
       geom_errorbar(aes(ymax=up,ymin=down)) +
       geom_line(aes(y=actual))

print(g)


load(paste0(datdir,'truebeta.RData'))
preds$cohort2 = preds$cohort^2
preds$cohort3 = preds$cohort^3
print(
  coefficients(lm(est~cohort+cohort2+cohort3,data=preds))
)

print(t.beta)


