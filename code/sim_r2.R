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
range(a.eff)
plot(tdat$a,a.eff)

draw_blocks = function(breakprob,autocorr,var){
    #breakprob is probability of a window break between 0 and 1
    #autocorr is an SD of effect sizes for a random walk
    #var is the variable to draw across rndomly
    
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
    var.beta=var.beta,var.eff=var.eff))

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


r2=0.03

evar=(1-r2)/r2
e=rnorm(nrow(tdat),0,sqrt(evar))

tdat$y = tdat$yhat + e
#tdat = dat
#tdat = dat %>% dplyr::select(-c,-a2,-p2,-c2)

print(
  summary(lm(y~a+p,data=tdat))
)

est = apcsamp(dat=tdat,dv='y',
              cores=4,method='ml',
              chains=4,samples=250)

summary(est)

effs = draw_effs(est)

#save(tdat,file=paste0(datdir,'nsim.RData'))
#pltdat=list()

#pull unique data -- should make this a loop
uns = list(
  c = unique(dat$c), a=unique(dat$a), p=unique(dat$p)
)

#need to hold means for margins 
dims=c('a','p','c')

for(d in dims){
  o.dims = dims[which(!dims==d)]
  o.means=apply(dat[,o.dims],2,mean); 
  xdat=data.frame(uns[[d]],o.means[1],o.means[2],
                  row.names=NULL); 
  colnames(xdat)=c(d,o.dims)
  xdat$a2 = xdat$a^2; xdat$p2 = xdat$p^2; xdat$c2 = xdat$c^2
  x = as.matrix(xdat)
  b = as.matrix(t.beta)[,colnames(x)]
  
  #test for matching
  if(!all(colnames(x) == colnames(b))){
    stop('Problem with matrix ordering.')}
  
  pltdat[[d]] = data.frame(
    id = uns[[d]],
    s1 = (x %*% b)[,1] 
  )
  
}

save(pltdat,file=paste0(datdir,'nsim_fits.RData'))

#run algorithm
#source('algorithm.R',echo=TRUE)
#unnecessary with "simulator.R"
