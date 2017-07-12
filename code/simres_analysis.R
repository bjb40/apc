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

simres = simres[extract('samp')>20]

print(unlist(extract('acc')))
extract('acc')

library(reshape2)
pvals=extract('pvals',as.df=FALSE,unlist=FALSE)
hist(unlist(pvals))

opv = do.call(c,lapply(pvals, function(l) l$omnibus))
besto = which(abs(opv-.5) == min(abs(opv-.5)))

rpv = do.call(rbind,lapply(pvals, function(l) range(l)))
rpv = cbind(rpv,rpv[,2]-rpv[,1]); colnames(rpv) = c('min.pval','max.pval','range.pval')
bestr = which(rpv[,3]==min(rpv[,3]))

pv = round(unlist(pvals),1)
sum(pv==.5)/length(pv)
sum(pv==.2)/length(pv)
avp=unlist(lapply(pvals,function(s) mean(unlist(s))))
best=which((abs(avp-0.5))==min(abs(avp-0.5)))
worst=which((abs(avp-0.5))==max(abs(avp-0.5)))

pv = melt(do.call(rbind,lapply(simres, function(s) s$pvals$a)))
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
  df$up = df$up - df$up[int]
  df$down=df$down-df$down[int]
  df$m_down=df$m_down-df$m_down[int]
  df$m_up = df$m_up-df$m_up[int]
  
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
  geom_point(aes(y=m_est),shape=1) +
  geom_ribbon(aes(ymin=m_down,ymax=m_up), alpha=0.25) +
  facet_grid(.~type,scales='free_x') +
  theme_classic()

return(plt)
}

#print(plt_fit(5))

#print(plt_fit(worst)); print(plt_fit(12)); print(plt_fit(13)); 
print(plt_fit(13)) #yuck --- print(plt_fit(12)) // has good omnibus pval --12 has bad!!
print(plt_fit(14)) ##14 and 15 look good! --- i wonder if rhat would help that... 
#23 has great example of best fit versus average!!!...


  