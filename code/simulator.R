#runs 5 simulations

#####check for and create simulations folder
source('config~.R')
if(!dir.exists(paste0(outdir,'simdata1/'))){
  dir.create(paste0(outdir,'simdata1/'))
}

####
#set up list of conditions for simulations

#simulate across blocks of years for APC

#random window blocks for P and C
yrblocks = c(3,5,10)
#simulation types for P and C
dim=c('p','c')
simtypes = c('ar1','random','randomwalk')
ac = c(0.2,0.5,0.8) #autocorrelations for ar1
r2 = 0.1 

simtable = expand.grid(yrblocks,yrblocks,
                       simtypes,simtypes,
                       ac,ac,r2)
names(simtable) = c('p.blocks','c.blocks',
                    'p.type','c.type',
                    'p.autocorr','c.autocorr','r2')
simtable = simtable %>%
  mutate(p.autocorr=ifelse(p.type !='ar1',NA,p.autocorr),
         c.autocorr=ifelse(c.type !='ar1',NA,c.autocorr))

simtable=unique(simtable)

simtable$nsim = as.numeric(NA)
simtable$ol = as.numeric(NA)
simtable$rrmse = as.numeric(NA)
simtable$srmse = as.numeric(NA) #ratio of calculated rmse to est rmse
simtable$p.catch = as.numeric(NA)
simtable$c.catch = as.numeric(NA)
#1:nrow(simtable)
for(draw in 1:5){
  simdims = simtable[draw,]
  cat('Simulation',drawsimulations,'...\n')
  #funciton is nested because sources below delete all in memory
  source('config~.R')
  source('sim_r2.R') #draws the random data / returns list
  simtable[draw,'rrmse'] = sqrt(evar)
  
  source('evalsim.R',echo=TRUE)
  simtable[draw,'ol'] = sum(mm$ol)/nrow(mm)
  simtable[draw,'srmse'] = mean(effs$fit$sigma)
  simtable[draw,'p.catch'] = 1-p.miss
  simtable[draw,'c.catch'] = 1-c.miss
  
  save(list = c('effs','simulated'),
       file=paste0(outdir,'simdata1/sim',drawsimulations,'.RData'))
  
  if(draw==1){append=FALSE} else{append=TRUE}
  
  write.table(simtable[draw,],
    file=paste0(outdir,'simdata1/simtable.csv'),
    append=append,
    col.names=!append,
    row.names=FALSE,
    sep=',')
  
}

