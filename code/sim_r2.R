#Bryce Bartlett
#simulates data for use in apc models

#clear cache
#rm(list=ls())
source('config~.R')

###############
#DRAW DATA
##############

#####
#SET DIMENSIONS FOR A and P (c wil be derived)

#preliminaries
library(dplyr)
library(apcwin)
library(parallel)

###
#draw random set of ages over a 
#certain number of periods (i.e.) 
#survey years (30) based on GSS 
p=seq(1972,2014,by=2)
#loop over p to create a draw of a,p,and c
tdat = lapply(p, function(y) {
  a=round(rnorm(2500,mean=46,sd=17.5)) # gss mean
  a = a[a>17 & a<79]
  df = data.frame(a=a,p=rep(y,length(a)))
  df$c=df$p-df$a
  return(df)
})
#bind
tdat = do.call(rbind,tdat)

n=nrow(dat)

####
#draw random effects values for various functions
tdat$a2 = (tdat$a-mean(tdat$a))^2
tdat$ac = tdat$a-mean(tdat$a)
#age effects = quadratic / centered on mean age
age = runif(1,-1,1)/10; age2=runif(1,-1,1)/10

a.eff = as.matrix(tdat[,c('ac','a2')]) %*% 
  matrix(c(age,age2),2,1)
#range(a.eff)
#plot(tdat$a,a.eff)

draw_blocks = function(breakprob,autocorr,var){
    #breakprob is probability of a window break between 0 and 1
    #autocorr is an SD of effect sizes for a random walk
    #var is the variable to draw across rndomly
  
  
    ####autocorr ideas : https://stats.stackexchange.com/questions/29239/creating-auto-correlated-random-values-in-r
    ####also https://stats.stackexchange.com/questions/178014/generating-an-ar1-time-series-with-a-specific-autocorrelation-in-r
    #### this works well
    "
    x=arima.sim(model=list(ar=0.5),n=100)
    acf(x)
    cor(cbind(x,stats::lag(x)),use='pairwise.complete.obs')
    "
  
    #returns list of (1) betas and (2) calculated effects for data in var
  
  #period and cohort effects have probability of "jump" 
  #and "autocorrelation"
  #breakprob = 0.3
  #autocorr = diff(range(a.eff))/5
  all.var = unique(var)
  #pull random breaks
  breaks.var = c(min(var)-1,
               all.var[as.logical(
                 rbinom(length(all.var),1,breakprob)
                 )])
  #make sure it covers the whole place
  if(max(breaks.var)<max(var)){
    breaks.var = c(breaks.var,max(var))
  }
  #draw random disturbance
  var.beta = rnorm(length(breaks.var)-1,sd=autocorr)
  #random walk from beginning
  var.beta = cumsum(var.beta)
  
  ##calculate effects
  #window is a method for apcwin that creates a factor
  vv=data.frame(v=window(var,breaks=breaks.var))
  var.eff = model.matrix(~v-1,data=vv) %*% 
    matrix(var.beta,length(var.beta),1)
  
  #rename
  names(var.beta) = levels(vv$v)
  
  return(list(
    var.beta=var.beta,var.eff=var.eff,breaks=breaks.var))

}

pdraw = draw_blocks(breakprob=0.3,
                    autocorr=diff(range(a.eff))/5, #1/5 of age effect
                    var=tdat$p)

period=pdraw$var.beta
plot(tdat$p,pdraw$var.eff)

cdraw = draw_blocks(breakprob=0.3,
                    autocorr=diff(range(a.eff))/5,
                    var=tdat$c)

cohort=cdraw$var.beta
plot(tdat$c,cdraw$var.eff)

tdat$yhat = a.eff + pdraw$var.eff + cdraw$var.eff


r2=0.1

#########
#r2 = 1- var(e)/var(y)
#var(e) = (1-r2)*var(y)
  #when var(y) = var(y_hat), var(e) = 0, and R^2 = 1 
  #so var(y) is the same as var(y_hat) + var(e), that means
#var(e) = (1-r2)*(var(y_hat) + var(e)) --- below, getting rid of variance
#e = y_hat + e - r2*y_hat - r2*e
#0 = y_hat -r2*y_hat - r2*e
#r2*e = y_hat +r2*y_hat
#e = (y_hat -r2*y_hat)/r2 

yhat = var(tdat$yhat)
evar=(yhat - r2*yhat)/r2
e=rnorm(nrow(tdat),0,sd=sqrt(evar))

tdat$y = tdat$yhat + e

tdat$pf = window(tdat$p,breaks=pdraw$breaks)
tdat$cf = window(tdat$c,breaks=cdraw$breaks)

#model
tt = lm(y~a+a2+pf+cf,data=tdat)
print(summary(tt)$r.squared)

####
#list for saving
simdat = tdat %>% dplyr::select(a,p,c,y)
names(pdraw$var.beta) = paste0('p',names(pdraw$var.beta))
names(cdraw$var.beta) = paste0('c',names(cdraw$var.beta))
age = c(age,age2); names(age) = c('cage','cage2')
effects = c(age,pdraw$var.beta,cdraw$var.beta)

##################33
#list of things to keep

simulated = list(
  simdat=simdat,
  effects=effects,
  p.breaks = pdraw$breaks,
  c.breaks = cdraw$breaks
)



