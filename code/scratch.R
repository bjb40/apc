########
#for new dgp


est = apcsamp(dat=tdat,dv='y',
              cores=4,method='ml',
              chains=4,samples=25)

summary(est)

use.marginal=TRUE

effs = draw_effs(est,tol=0.1,marginal=use.marginal)



###
#real data

gm = grand.means=apply(tdat[,c('a','p','c')],2,mean)

#make p and c into mean level
gm['p'] = window(gm['p'],breaks=pdraw$breaks)
gm['c'] = window(gm['c'],breaks=cdraw$breaks)


#pool real effects -- add grand mean of apc
r.effs = list(a=a.eff, 
              p=pdraw$var.eff,
              c=cdraw$var.eff)

if(use.marginal){
  r.effs = list(a=a.eff + period[gm['p']] + cohort[gm['c']], 
                p=pdraw$var.eff + gm['a']*age + gm['a']^2*age2 + cohort[gm['c']],
                c=cdraw$var.eff + gm['a']*age + gm['a']^2*age2 + period[gm['p']])
}

#####
#calculate effects
xhat=data.frame(
  a=unique(tdat$a)
)
xhat$a = xhat$a^2


#pull unique values for plotting
r.effs=lapply(names(r.effs),function(x){
  df = data.frame(eff=r.effs[[x]],val=tdat[,x])
  df$dim=x
  return(unique(df))
})  

r.effs = do.call(rbind,r.effs)

pp = plot(effs)

pp +
  geom_line(data=r.effs,
            aes(x=val,y=eff),linetype=2) + 
  theme_classic()

#ggplot(r.effs,aes(x=val,y=eff)) + 
#  geom_line(linetype=2) + 
#  theme_classic() +
#  facet_wrap(~type,scales='free_x')

mm = merge(r.effs %>% rename(x=val),pp$data,by=c('dim','x'))

mm$delta = mm$eff - mm$Fit

mm %>% 
  group_by(dim) %>%
  summarize(eff=mean(eff),
          fit=mean(Fit),
          delta=mean(delta))

mm=mm %>% group_by(dim) %>%
  mutate(fit.dm = Fit - mean(Fit))



#model is identifiable!! (can check directly from betas..., but using scopedummy and predict)
tt = lm(y~a+a2+pf+cf,data=tdat)
pred=list()

###a effects
xh= data.frame(a=min(tdat$a):max(tdat$a)); xh$a2 = xh$a^2
xh$pf = window(mean(tdat$p),breaks=pdraw$breaks)
xh$cf = window(mean(tdat$c),breaks=cdraw$breaks)

pred[['a']] = data.frame(
  x=xh$a,
  m.eff = predict(tt,newdata=xh))
pred[['a']]$dim='a'

###p effects
xh.c= data.frame(cf=scopedummy(tdat$cf))
xh.c$a = mean(tdat$a); xh$a2=mean(tdat$a2)
xh.c$pf = window(mean(tdat$p),breaks=pdraw$breaks)


##c effects
pred[['c']] = data.frame(
  x=min(tdat$c):max(tdat$c),
  m.eff = predict(tt,newdata=xh))



tst=merge(mm %>% filter(dim=='c'),pred,by='x')
ggplot(tst,aes(x=x,y=eff)) + 
  geom_line(color='red',linetype=2) + 
  geom_line(aes(x=x,y=m.eff),color='blue',linetype=2) +
  geom_line(aes(x=x,y=Fit)) +
  geom_ribbon(aes(ymax=ul,ymin=ll),alpha=0.25) + 
  theme_classic()



ggplot(mm,aes(x=x,y=Fit)) +
  geom_line(color='blue') + 
  geom_line(aes(x=x,y=fit.dm,color='red')) + 
  geom_line(aes(x=x,y=eff),linetype=2) +
  facet_grid(~dim,scales='free')


grand.means = data.frame(
  a = window(mean(tdat$a),breaks=breaks$a[[s]]),
  p = window(mean(tdat$p),breaks=breaks$p[[s]]),
  c = window(mean(tdat$c),breaks=breaks$c[[s]])
)

