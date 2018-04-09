######
#bind simtable with simulator number
#simulator.R did not include a reference to simulation number
#this corrects the problem!

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

simtable$simnumber = 1:nrow(simtable)

ss = read.csv('H:/projects/apc/output/simdata1/simtable.csv')

mrg = merge(simtable,ss,all=TRUE)


#add additional columns for HAPC tests
mrg$hol = NA #overlap of CI for simulated HAPC
mrg$hrmse = NA #sigma for HAPC
mrg$efficiency = NA #this is efficiency of our approach rel to HAPC


write.csv(mrg,'H:/projects/apc/output/simdata1/fullsimtable.csv',
          row.names=FALSE)
