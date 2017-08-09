#analyze simulation data

#########
#NOTE THAT YOU HAVE TO FIX THE PROBLEM ---- THE FIRST ONES (1-10) DON'T ACTUALLY CHANGE TBETA wen del
#they are duplicates, also
####

source('config~.R')

#extraction function --- using names
extract = function(name,as.df=TRUE,unlist=FALSE){
  #input is a name of the simres object
  #output is a numbered list (by simulaiont)
  r=lapply(simres, function(s) s[[name]])
  if(as.df){r=do.call(rbind,r)}else if(unlist){r=unlist(r)}
  return(r)
}


load(paste0(outdir,'simres.RData'))

#use first 6 sufficient simulations
simres=simres[extract('samp')>=1000 & extract('bound')<.25]
simres=simres[1:6]

print(unlist(extract('acc')))
extract('acc')

library(reshape2)
pvals=extract('pvals',as.df=FALSE,unlist=FALSE)
hist(unlist(pvals))

opv = do.call(base::c,lapply(pvals, function(l) l$omnibus))
besto = which(abs(opv-.5) == min(abs(opv-.5)))
worsto = which(abs(opv-.5) == max(abs(opv-.5)))

rpv = do.call(rbind,lapply(pvals, function(l) range(l)))
rpv = cbind(rpv,rpv[,2]-rpv[,1]); colnames(rpv) = c('min.pval','max.pval','range.pval')
bestr = which(rpv[,3]==min(rpv[,3]))

pv = round(unlist(pvals),1)
sum(pv==.5)/length(pv)
sum(pv==.2)/length(pv)
avp=unlist(lapply(pvals,function(s) mean(unlist(s))))
best=which((abs(avp-0.5))==min(abs(avp-0.5)))
worst=which((abs(avp-0.5))==max(abs(avp-0.5)))

pv = melt(do.call(rbind,lapply(simres, function(s) s$pvals$p)))
#box plot of p-values
hg = ggplot(pv,aes(x=factor(Var2),y=value)) + 
  geom_boxplot() +
  geom_point(aes(y=0.5,x=factor(Var2)), shape=8,color='grey',size=.25) +
  theme_classic()

print(hg)


###########

#have to fix recenter.... b/c the error bars shouldn't change by avg, but relative to
#the mean... would be better if I had a diff.....
recenter=function(df){
  
  delts = df %>% 
    mutate(up = up-est, 
           down=down-est, 
           m_down = m_down-m_est,
           m_up = m_up-m_est) %>%
    dplyr::select(up,down,m_up,m_down,id)
  
  rownames(delts) = rownames(df)
  
  mn = df %>%
    summarize_all(mean)

  tst = as.data.frame(sweep(as.matrix(df),2,as.matrix(mn),FUN='-')) %>% 
    dplyr::select(actual,est,m_est,id)
  tst$id=tst$id + mn$id

  tst = merge(tst,delts,by='id')
  
  tst = tst %>% 
    mutate(up=up+est,
           down=down+est,
           m_down=m_down+m_est,
           m_up = m_up+m_est)
  
  return(tst)
}


plt_fit = function(simnum,recenter=FALSE,legend=FALSE){
  #input simulaiton number; output is graph of effects
pred = simres[[simnum]]$preds
if(recenter){pred=lapply(pred,recenter)} #need fixed b/c first is not always continuous
df=do.call(rbind,pred)
df$type=substr(rownames(df),1,1)

if(legend){
  ls = c('1' = 'Actual', '2' = 'Estimated')
  
df$lty1 = as.factor(ls[1])
df$lty2 = as.factor(ls[2])
plt=ggplot(df,aes(x=id,y=actual,linetype=lty1)) + 
  labs(linetype='') +
  geom_line(aes(y=m_est,linetype=lty2))
  
} else{
  plt=ggplot(df,aes(x=id,y=actual)) +
    geom_line(aes(y=m_est))
  
}

plt = plt + geom_line() +
  #geom_point(aes(y=est), alpha=0.25) +
  #geom_errorbar(aes(ymin=down,ymax=up), alpha=0.5) +
  geom_ribbon(aes(ymin=m_down,ymax=m_up), alpha=0.25) +
  facet_grid(.~type,scales='free_x') +
  theme_minimal() +
  ylab('Effect') +
  xlab('APC Value') 


return(plt)
}

#print(plt_fit(5))

#print(plt_fit(worst)); print(plt_fit(12)); print(plt_fit(13)); 
#print(plt_fit(13)) #yuck --- print(plt_fit(12)) // has good omnibus pval --12 has bad!!
#print(plt_fit(14)) ##14 and 15 look good! --- i wonder if rhat would help that... 
#21 & 15 has great example of best fit versus average!!!...

label_add = labs(title='Results for Bayesian Model Averaging Window Constraint APC Models.',
       caption='n=1,000\nAverage of 1,000 Different Block Models Sampled Using MCMC.')

#print(plt_fit(21))
print(plt_fit(worsto, recenter=TRUE,legend=TRUE))
ggsave(paste0(imdir,'worst_sim.pdf'))
print(plt_fit(besto, recenter=TRUE, legend=TRUE))
ggsave(paste0(imdir,'best_sim.pdf'))

###twitter

print(plt_fit(besto,recenter=TRUE,legend=TRUE)+label_add)
ggsave(paste0(imdir,'twitter_sim.png'),
       height=8,width=12,units='in',dpi=150)

######
#table of overlaps...
######


ol.preds = function(preds,recenter=FALSE){
  #calculates overlap of preds --
  #returns proportion of actual within 95% intervals
  if(recenter){preds=lapply(preds,recenter)} #need fixed b/c first is not always continuous
  preds = do.call(rbind,preds)
  ol = preds$m_up > preds$actual & preds$m_down < preds$actual
  return(sum(ol)/length(ol))
}

ol.unadj = unlist(lapply(simres,function(x) ol.preds(x$preds)))
ol.adj = unlist(lapply(simres,function(x) ol.preds(x$preds,recenter = TRUE)))

opvd = abs(.5-opv)
avpd = abs(.5-avp)
rpvd = rpv[,'range.pval']

fdat = data.frame(ol.adj,ol.unadj,opvd,avpd,rpvd,
                  samp=extract('samp'),bound=extract('bound'),
                  r2=extract('r2')[,'Median'],acc=extract('acc'))

r = fdat %>% 
  dplyr::select(opvd,avpd,rpvd) %>% 
  mutate(rpvd=1/rpvd)

r=apply(r,2,rank)

fdat$rank = rowMeans(r[,1:2])
  
  
print(round(cor(fdat),3))

fdatm = melt(fdat %>% dplyr::select(-samp,-ol.unadj),id='ol.adj')

ggplot(fdat,aes(y=ol.adj,x=opvd,color=acc)) +
  geom_point() +
  theme_classic()

plts2=lapply(1:6,plt_fit,recenter=TRUE)
plts2 = lapply(plts2,function(x) x+theme_void())

library(gridExtra)

  grid.arrange(grobs=plts2,nrow=3)

##output table
library(knitr)

tb = extract('tbeta')
