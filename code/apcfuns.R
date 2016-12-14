#Dev R 3.3.0 "Supposedly Educational"
#general functions for apc project
#
#Bryce Bartlett


library(ggplot2)

#make an S3 object: 
#http://www.cyclismo.org/tutorial/R/s3Classes.html#memory-management

window = function(var,winlength,breaks){
  #this function makes windows from a continuous meaasures
  #
  #Input: var = numeric vector
  #       winlength = number identifying width of windows
  #       breaks = vector of breaks (like cut in base)
  #       can only provide 1 winlenth or breaks
  
  #check and throw errors
  if(missing(winlength)){
    winlength=NA
  }
  
  if(missing(breaks)){
    nobreak=TRUE
  } else{nobreak=FALSE}
  
  if(nobreak & is.na(winlength)){
    stop("
Must specify either a window length (winlength)
or vector of breaks (breaks -- works like base cut)")
  }
  
  #consider warning for small cell values
  #consider error for out-of-range window or break sets
  
  r=range(var)
  
if(!is.na(winlength)){
  w=winlength
  #full possibility
  vec=r[1]:r[2]
  #vector of windows
  win=1:w
  breaks = r[1]-1
  for(i in 0:(ceiling(length(vec)/w)-1)){
    p.break=vec[max(win+w*i)]
    if(is.na(p.break)){p.break=r[2]}
    breaks=append(breaks,p.break)
  }
}#end winlength 

  c=cut(var,breaks=breaks)
  attr(c,'range') = r
  attr(c,'breaks') = breaks
  #inherit the name???
  
  class(c) = append(class(c),'window')
  
  return(c)
  
}

scopedummy=function(w){
  #this is a method for window object that transforms it
  #into a factor covering all possible continuous values
  #the (row) names are the values
  #i.e. 1,2,3,4,5,6 with two windows output two factors
  # 1,1,1,2,2,2
  
  UseMethod('scopedummy',w)

}

scopedummy.default = function(w) {
  cat('\nError: Must use window object.\n\n')
}

scopedummy.window = function(w) {
  r=attr(w,'range')
  span=r[1]:r[2]
  f=cut(span,breaks=attr(w,'breaks'))
  return(f)
}

range=function(w){
  UseMethod('range',w)
}

range.window = function(w){
  #returns range object (from underlying continuous)
  return(attr(w,'range'))
}

#ols funciton for estimating an apc model
apc_lm = function(formula,data,age,per,coh,windows) {
  #runs a linear model on an apc problem (need to id window)
  #
  # Inputs
  #  formula=r formula object
  #  data = dataframe
  #  age,per,coh = named columns of dataframe for age, period, cohort
  #  need only specify 2... If all 3 are specified, model checks for
  #  linear dependence
  #  windows = list identifying window constraints across apc
  #            if empty, will pick random windows min=3,max=max(variable))
  #
  # Output
  #   list including 
  #   (1) effects =  ols estimates from window model
  #   (2) smallblock = block APC estimates
  #   (3) linear = linearlized block estimates (cubic)

  #@@@@
  #input checks
  no.age=no.period=no.cohort=FALSE
  
  if(missing(age)){
    no.age=TRUE
    age='age'
  } 
  if(missing(per)){
    no.period=TRUE
    per='period'
  }
  if(missing(coh)){
    no.cohort=TRUE
    coh='cohort'
  }
  
  if(sum(no.period,no.cohort,no.age)>1){
    stop('Must specify 2 of 3 required estimates: age, period, or cohort.')
  } else if(no.age){
      data[,age]=data[,per]-data[,coh]
  } else if(no.period){
      data[,per]=data[,age]+data[,coh]
  } else if(no.cohort){
      data[,coh]=data[,per]-data[,age]
  }
  
  
  
  #@@@@
  #window check
  if(missing(windows)){
    windows=list(age=0,period=0,cohort=0)
    
    #id maximum window for constraint(min is 3)
    max_w=function(var){
      r=range(var)
      return(r[2]-r[1])
      }
    windows$age=round(runif(1,3,max_w(data[,age])))
    windows$period=round(runif(1,3,max_w(data[,per])))
    windows$cohort=round(runif(1,3,max_w(data[,coh])))
    
  }
  
    
  #@@@@
  #build model matrix from window constraints
  wins=list(
    a=window(data[,age],winlength=windows$age),
    p=window(data[,per],winlength=windows$period),
    c=window(data[,coh],winlength=windows$cohort)
  )
  

  ndat=data[,!colnames(data) %in% c(age,per,coh)]
  ndat=cbind(ndat,wins$a,wins$p,wins$c)
  colnames(ndat)[-1:(3-ncol(ndat))] = c(age,per,coh)
  
  blockdat=lapply(wins,scopedummy)
  
  #@@@@
  #estimate OLS model
  
  results=lm(formula,data=ndat)
  
  #@@@@
  #prepare small block estimates
  
  #x.per=blockdat$per
  
  #@@@@
  #prepare linear estimates(cubic)
  
  #@@@@
  #return
  
  return(
    list(
      newdat=ndat,
      blockdat=blockdat,
      results=results
    )
  )
  
}

plt = function(ols,varnames){
  #this function takes a,p,c measures
  #and plots them in a panel of 3
  #dependency -- ggplot2 (can remove for package)
  #
  #Input: 
  #     ols=lm object
  #     varnames=charater vector (1-3) identifying
  #     variablenames for age period and cohort (in that order)
  #     This will work for prefixes...
  #
  #Output: list of ggplot objects
  
  #select coefficients from 
  b = coefficients(ols)
  
  return(b)
  
}

