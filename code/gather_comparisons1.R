###
#this script gathers comparisons from tst_comparisons
source('config~.R')
#dependencies
library(apcwin) #custom package in dev
library(lme4) #random effects
library(merTools) #random effects display and sampling
library(parallel) #to speed up sampling
library(dplyr) #data manipulation
library(reshape2) # data manipulation
library(ggplot2) #plotting

simdir = paste0(outdir,'simdata1/')
resfile = paste0(simdir,'comp.csv')
simtable = read.csv(paste0(simdir,'fullsimtable.csv')) %>%
  arrange(simnumber)

sims = simtable$simnumber

for(s in 1:length(sims)){
  #load data
  load(paste0(simdir,'sim',sims[s],'.RData'))
  #calculate revised comparisons
  source('tst_comparisons.R')
  #merge into "revtable"
  overview['simnumber'] = sims[s]
  
  #save comparisons
  append=file.exists(resfile) 
  write.table(t(overview),
              file=resfile,
              append=append,
              col.names=!append,
              row.names=FALSE,
              sep=',')
  
  #update RData with new objects
  truemodel = tt; rm(tt)
  
  save(list=c('effs','hapc','truemodel','summdat','simulated'),
    file=paste0(simdir,'sim',sims[s],'.RData'))
  
}

  
 