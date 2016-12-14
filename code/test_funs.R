#simulates data for 4 scenarios from Luo & Hodges SMR 2015; p. 14
#
#Bryce Bartlett

#clear cache
rm(list=ls())
source('config~.R')

library(microbenchmark)

#load test data
load(paste0(datdir,'testdat.RData'))

testwindow=FALSE

if(testwindow){
######
#WINDOW OBJECT TESTS
######

a=tdat$a

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
a = tdat$a
w = window(a,winlength=5)

#should throw error:
scopedummy(a)

invert.window=scopedummy(w)

print(
  table(invert.window)
)

######
#ols linear model tests
######




