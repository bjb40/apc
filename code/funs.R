#Dev R 3.3.0 "Supposedly Educational"
#general functions for apc project
#
#Bryce Bartlett


library(ggplot2)

setwindow = function(var,w){
  #this funciton makes windows from a continuous meaasures
  #
  #Input: var = numeric vector
  #       w = number identifying width of windows
  
  return(cut_interval(var,length=w))
  
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

#test
testb = e.beta[[1]]
plt(testb,varname=c('a','p','c'))
