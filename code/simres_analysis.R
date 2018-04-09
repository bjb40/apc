#analyze simulation data

#########
#NOTE THAT YOU HAVE TO FIX THE PROBLEM ---- THE FIRST ONES (1-10) DON'T ACTUALLY CHANGE TBETA wen del
#they are duplicates, also
####

source('config~.R')
library(scales)

#load simulation results

simtable = read.csv(
  paste0(outdir,'simdata1/fullsimtable.csv'))


ggplot(simtable,
       aes(x=c.autocorr,y=p.autocorr,color=efficiency,
           shape=factor(efficiency>1))) +
  scale_color_gradient2(high=muted('red'),mid='white',low=muted('blue'),
                        midpoint=1) +
  geom_jitter(width=0.05,height=0.05) +
  theme_classic()


ggplot(simtable,
       aes(x=c.blocks,y=p.blocks,color=efficiency,
           shape=factor(efficiency>1))) +
  scale_color_gradient2(high=muted('red'),mid='white',low=muted('blue'),
                        midpoint=1) +
  geom_jitter(width=0.5,height=0.5) +
  theme_classic()

tst = melt(simtable %>% 
  ungroup %>% 
    dplyr::select(efficiency,p.type,c.type),id='efficiency')

ggplot(tst,
       aes(x=value,y=efficiency,color=variable)) +
  geom_boxplot() +
  theme_classic()

hist(simtable$efficiency) #not quite normal...

tst = simtable[,c(1:4,11,13:14,17)]
summary(lm(efficiency~.,data =tst))
summary(lm(efficiency~p.blocks*p.type + c.blocks*c.type,data =tst))
summary(lm(efficiency~factor(p.autocorr)*p.blocks + 
             factor(c.autocorr)*c.blocks,data=simtable))

#summary(lm(efficiency~p.autocorr*p.blocks + c.autocorr*c.blocks,data=simtable))

