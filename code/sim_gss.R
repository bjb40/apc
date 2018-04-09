#Bryce Bartlett
#simulates data for use in apc models

#################333
#note that this refers to variables built in sim_r2!!

#clear cache
#rm(list=ls())
source('config~.R')

####
#these need given
#simdims = list(p.blocks = 5, c.blocks = 5,
#               p.type='random',c.type='random',
#               p.autocorr=NA,c.autocorr=NA
#               )

useGSS=simdims$useGSS

###############
#DRAW DATA
##############

#####
#SET DIMENSIONS FOR A and P (c wil be derived)

#preliminaries
library(dplyr)
library(apcwin)
library(parallel)


####
#load a,p,c from GSS cross section
dat.f = 'H:/projects/proposal/r_study/output/private~/cxgss.RData'


if(file.exists(dat.f)){
  load(dat.f)
  #save object for runing on the server
  save(cleandat,file=paste0(outdir,'~dat_parsol.RData')) #~ keeps it from pushing to git
} else{
  load(paste0(outdir,'~dat_parsol.RData'))
}

gssdat = cleandat %>%
  rename(a=age,
         p=year,
         c=birthyear) %>%
  dplyr::select(a,p,c) 
rm(cleandat)
gssdat = gssdat[complete.cases(gssdat),]

###
#draw random set of ages over a 
#certain number of periods (i.e.) 
#survey years (30) based on GSS 

if(!useGSS){
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
} else if(useGSS){
  tdat = gssdat
}


#n=nrow(tdat)

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

#need to change the variance

ar1 = function(len,autocorr,v){
  #this simulates a stationary AR1 process with correlation autocor
  #len is length
  #autocor is autocorrelation
  #v is variance: see https://stats.stackexchange.com/questions/175289/how-to-set-the-starting-value-of-simulated-arima
  #for calculation of v
  
  var.beta = arima.sim(model=list(ar=autocorr),
                       n=len)
  
  #vratio = v/var.beta[1]
  var.beta=var.beta*as.vector(sqrt(v))
  
  return(var.beta)
}

######
#unit tests 
#to confirm stationary conttruction with wider variance
"tst1 = ar1(100,0.5,1)
  mean(tst1); var(tst1)
  print(acf(tst1,plot=FALSE)[1])
tst10 = ar1(100,0.5,10)
  mean(tst10); var(tst10)
  print(acf(tst10,plot=FALSE)[1])
"

randomwalk = function(len,v){
  #this simulates a random walk based on the last observation
  #len is length and v is variance
  #NOTE that random walks have fairly high autocorrelation by nature
  #https://stats.stackexchange.com/questions/181167/what-is-the-autocorrelation-for-a-random-walk
  var.beta = cumsum(rnorm(len,sd=sqrt(v)))
  return(var.beta)
}

random = function(len,beg,end){
  #this just draws from a uniform random distribution beginning and end
  #this is a function for uniformity with the above
  return(runif(len,beg,end))
}

draw_blocks = function(breakyears=5,
                       type='random',
                       autocorr=0,
                       hetvar=FALSE,
                       var){
    #breakyears is average length of break years (more or less)
    #type is one of 
    #  random-walk, ar1, or random 
      #in all cases, effect sizes are roughly the same range as age effects
    #autcor is autocorrelation number between 0 and 1 for stationary process
    #hetvar is the size of a heterscedastic variance; 0 is no heterscedastic variance
    #a number indicates an additive or subtractive variance
    #var is the variable to draw across rndomly
  
  
    ####autocorr ideas : https://stats.stackexchange.com/questions/29239/creating-auto-correlated-random-values-in-r
    ####also https://stats.stackexchange.com/questions/178014/generating-an-ar1-time-series-with-a-specific-autocorrelation-in-r
    #### this works well

  #period and cohort effects have probability of "jump" 
  #and "autocorrelation"
  #breakprob = 0.3
  #autocorr = diff(range(a.eff))/5
  yrs = min(var):max(var)
  breakprob = length(yrs)/breakyears/length(yrs)
  all.var = unique(var)
  #pull random breaks
  breaks.var = c(min(var)-1,
                 yrs[as.logical(
                 rbinom(length(yrs),1,breakprob)
                 )])
  #make sure it covers the whole place
  if(max(breaks.var)<max(var)){
    breaks.var = c(breaks.var,max(var))
  }
  
  #draw random disturbance---these are one of several choices, id'd and discussed in
  #functions above
  ###########
  if(type=='random'){
    var.beta = random(len=length(breaks.var)-1,
                      beg=min(a.eff),
                      end=max(a.eff))
  } else if(type=='randomwalk'){
    var.beta = randomwalk(len=length(breaks.var)-1,
                          v=var(a.eff)/2)
  } else if(type=='ar1'){
    var.beta = ar1(len=length(breaks.var)-1,
                   autocorr=autocorr,
                   v=var(a.eff))
  }
  
  ##calculate effects
  #window is a method for apcwin that creates a factor
  vv=data.frame(v=window(var,breaks=breaks.var))
  var.eff = model.matrix(~v-1,data=vv) %*% 
    matrix(var.beta,length(var.beta),1)
  
  if(hetvar){
    var.sd = randomwalk(len=length(breaks.var)-1,
                        v=(var(var.eff)/5))
    var.sd=model.matrix(~v-1,data=vv) %*%
      matrix(sqrt(abs(var.sd)),length(var.sd),1)
  
    var.eff = rnorm(n=length(var.eff),mean=var.eff,sd=var.sd)
  }
        
  #rename
  names(var.beta) = levels(vv$v)
  
  return(list(
    var.beta=var.beta,var.eff=var.eff,breaks=breaks.var))

}


###################
#draw random effects from pre-specified list simlist
#in simulator.R

pdraw = draw_blocks(breakyears=simdims$p.blocks,
                    type=simdims$p.type,
                    autocorr=simdims$p.autocorr, #1/5 of age effect
                    var=tdat$p)

period=pdraw$var.beta
#plot(tdat$p,pdraw$var.eff)

cdraw = draw_blocks(breakyears=simdims$c.blocks,
                    type=simdims$c.type,
                    autocorr=simdims$c.autocorr,
                    var=tdat$c)

cohort=cdraw$var.beta
#plot(tdat$c,cdraw$var.eff)

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



