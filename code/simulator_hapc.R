#runs simulations to test sensitivity relative to HAPC under 2 conditons
#1) low cell values
#2) non-stationary, non-uniform variance time series

#####check for and create simulations folder
source('config~.R')
if(!dir.exists(paste0(outdir,'simdata_hapc/'))){
  dir.create(paste0(outdir,'simdata_hapc/'))
}

####
#set up list of conditions for simulations

#simulate across blocks of years for APC

#random window blocks for P and C
yrblocks = c(3,5,10)
#simulation types for P and C
dim=c('p','c')
simtypes = c('ar1','random','randomwalk')
hetvar = TRUE
useGSS = TRUE
ac = 0.5 #autocorrelations for ar1
r2 = 0.1 

simtable = expand.grid(yrblocks,yrblocks,
                       simtypes,simtypes,
                       hetvar,useGSS,
                       ac,ac,r2)
names(simtable) = c('p.blocks','c.blocks',
                    'p.type','c.type',
                    'hetvar','useGSS',
                    'p.autocorr','c.autocorr','r2')
simtable = simtable %>%
  mutate(p.autocorr=ifelse(p.type !='ar1',NA,p.autocorr),
         c.autocorr=ifelse(c.type !='ar1',NA,c.autocorr))

simtable=unique(simtable)

simtable$simnumber = 1:nrow(simtable)
simtable$ol = as.numeric(NA)
simtable$rrmse = as.numeric(NA)
simtable$srmse = as.numeric(NA) #ratio of calculated rmse to est rmse
simtable$p.catch = as.numeric(NA)
simtable$c.catch = as.numeric(NA)

for(draw in simtable$simnumber){
  simdims = simtable[simtable$simnumber==draw,]
  cat('Simulation',draw,'...\n')
  #funciton is nested because sources below delete all in memory
  source('config~.R')
  source('sim_gss.R') #draws the random data / returns list
  simtable[draw,'rrmse'] = sqrt(evar)
  
  source('evalsim.R',echo=TRUE)
  simtable[draw,'ol'] = sum(mm$ol)/nrow(mm)
  simtable[draw,'srmse'] = mean(effs$fit$sigma)
  simtable[draw,'p.catch'] = 1-p.miss
  simtable[draw,'c.catch'] = 1-c.miss
  
  save(list = c('effs','simulated'),
       file=paste0(outdir,'simdata_hapc/sim',draw,'.RData'))
  
  if(draw==1){append=FALSE} else{append=TRUE}
  
  write.table(simtable[draw,],
    file=paste0(outdir,'simdata_hapc/simtable.csv'),
    append=append,
    col.names=!append,
    row.names=FALSE,
    sep=',')
  
}

