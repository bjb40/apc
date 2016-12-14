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
    breaks=append(breaks,max(win+w*i))
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

#ols funciton for estimating an apc model
apc_lm = function(formula,data,a,p,c,windows) {
  #runs a linear model on an apc problem (need to id window)
  #
  # Inputs
  #  formula=r formula object
  #  data = dataframe
  #  a,p,c = named columns of dataframe for age, period, cohort
  #  need only specify 2... If all 3 are specified, model checks for
  #  linear dependence
  #  windows = list identifying window constraints across apc
  #
  # Output
  #   list including 
  #   (1) effects =  ols estimates from window model
  #   (2) smallblock = block APC estimates
  #   (3) linear = linearlized block estimates (cubic)
  
  
  
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

