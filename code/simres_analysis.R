#analyze simulation data

#########
#NOTE THAT YOU HAVE TO FIX THE PROBLEM ---- THE FIRST ONES DON'T ACTUALLY CHANGE TBETA wen del
####

source('config~.R')
load(paste0(outdir,'simres.RData'))

#extraction function --- using names
extract = function(name,as.df=TRUE,unlist=FALSE){
  #input is a name of the simres object
  #output is a numbered list (by simulaiont)
  r=lapply(simres, function(s) s[[name]])
  if(as.df){r=do.call(rbind,r)}else if(unlist){r=unlist(r)}
  return(r)
}

print(unlist(extract('acc')))
extract('acc')

library(reshape2)
pvals=extract('pvals',as.df=FALSE,unlist=FALSE)
hist(unlist(pvals))

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

recenter=function(df){
  int = df$id == min(df$id)
  df$est = df$est - df$est[int]
  df$actual = df$actual - df$actual[int]
  df$m_est = df$m_est - df$m_est[int]
  
  return(df)
}


plt_fit = function(simnum){
  #input simulaiton number; output is graph of effects
pred = simres[[simnum]]$preds
#pred=lapply(pred,recenter)
df=do.call(rbind,pred)
df$type=substr(rownames(df),1,1)

plt=ggplot(df,aes(x=id,y=actual)) +
  geom_line(lty=2) + 
  geom_point(aes(y=est)) +
  geom_ribbon(aes(ymin=down,ymax=up),alpha=0.25) +
  #geom_point(aes(y=m_est),shape=2) +
  #geom_ribbon(aes(ymin=m_down,ymax=m_up), alpha=0.25) +
  facet_grid(.~type,scales='free_x') +
  theme_classic()

return(plt)
}

print(plt_fit(best)); #print(plt_fit(12)); print(plt_fit(13))
print(plt_fit(worst))



  