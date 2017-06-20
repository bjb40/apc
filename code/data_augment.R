#data_augmentation
#bryce bartlett
#this uses a data augmentation strategy to "smudge" the edges of the apc
#and break the identification problem


#clear cache
rm(list=ls())
source('config~.R')


#prelim
dv='y1'
actual='s1'

#load test data
#load(paste0(datdir,'testdat.RData'))
#load(paste0(datdir,'nsim.RData'))
#load(paste0(datdir,'sim2_tbeta.RData'))
load(paste0(datdir,'testdat.RData')) #luo and hodges
tdat$c = tdat$p-tdat$a

#smudge using uniform distributions -.5 to .5 ---> assuming they are closer to
#one or the other place

##helper functions

########
#smudge functions
########

smudge =function(var){
  return(var + runif(length(var),-.5,.5))
}

smudge_dat = function(df,vars){
  #df is dataframe
  #vars is indexes or character vector

  #select appropriate columns (as vars)
  if(missing(vars)){
    vars=1:ncol(df)
  }
    
  
  if(is.character(vars)){
    vars = colnames(df) %in% vars
  }
  
  if(is.numeric(vars)){
    vars = 1:ncol(df) %in% vars
  }
  
  
  df = cbind(df[,!vars],apply(df[,vars],2,smudge))
  return(df)

}

########
#MI rubins rules
#all of these operate on a list of lm fit objects
########

point_est = function(fits){
  est = data.frame(lapply(fits,coef)); colnames(est) = 1:ncol(est)
  return(rowMeans(est))
}

var_est = function(fits){
  var = data.frame(lapply(fits,function(x) diag(vcov(x))))
  est = data.frame(lapply(fits,coef)); colnames(est) = 1:ncol(est)
  
  colnames(var)=1:ncol(var)
  
  within = rowMeans(var)
  between = rowSums((rowMeans(est) - est)^2)/(length(fits)-1)
  total = within + between + (between/length(fits))
  
  return(list(
    within=within,
    between=between,
    total=total)
    )
  
}


tvals_est = function(fits){
  est = point_est(fits)
  var = var_est(fits)
  m=length(fits)
  
  tstat = est/var$total
  #tstat = sqrt(var$total)*(0-est)
  df.ratio=var$within/((1+1/m)*var$between)
  df=(m-1)*sqrt(1+df.ratio)

  pv = 2*(1-pt(abs(tstat),df=df))
  return(pv)
}

####
#example
####

n.sim = 150

dats = lapply(1:30,function(x) smudge_dat(tdat,c('a','p','c')))

fits = lapply(dats,function(x)
  lm(y3~a+p+c+I(a^2)+I(p^2)+I(c^2)-1,data=x))

print('mi fit')
print(point_est(fits))
print(tvals_est(fits))

print('real beta')
#print(t.beta[c('a','p','c','a2','p2','c2')])
#luo and hodges: a=0.3,a2=-0.01,p=-0.04,p2=0.02,c=0.35,c2=-0.0015

#print(mean(unlist(lapply(fits,function(x) summary(x)$r.squared))))

###plot some stuff

