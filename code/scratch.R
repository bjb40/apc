########
#for new dgp


est = apcsamp(dat=tdat,dv='y',
              cores=4,method='ml',
              chains=4,samples=25)

summary(est)

use.marginal=TRUE

effs = draw_effs(est,tol=0.1,marginal=use.marginal)


#model is identifiable!! (can check directly from betas..., but using scopedummy and predict)
tdat$a2 = tdat$a^2
tt = lm(y~a+a2+pf+cf,data=tdat)
true.b =list()

###a effects
xh= data.frame(a=min(tdat$a):max(tdat$a)); xh$a2 = xh$a^2
xh$pf = window(mean(tdat$p),breaks=pdraw$breaks)
xh$cf = window(mean(tdat$c),breaks=cdraw$breaks)

true.b[['a']] = data.frame(
  x=xh$a,
  m.eff = predict(tt,newdata=xh))
true.b[['a']]$dim='a'

###p effects
xh= data.frame(pf=scopedummy(tdat$pf))
xh$a = mean(tdat$a); xh$a2=mean(tdat$a2)
xh$cf = window(mean(tdat$c),breaks=cdraw$breaks)

true.b[['p']] = data.frame(
  x=min(tdat$p):max(tdat$p),
  m.eff = predict(tt,newdata=xh))
true.b[['p']]$dim='p'

###c effects
xh= data.frame(cf=scopedummy(tdat$cf))
xh$a = mean(tdat$a); xh$a2=mean(tdat$a2)
xh$pf = window(mean(tdat$p),breaks=pdraw$breaks)

true.b[['c']] = data.frame(
  x=min(tdat$c):max(tdat$c),
  m.eff = predict(tt,newdata=xh))
true.b[['c']]$dim='c'

#####
#plot

r.effs = do.call(rbind,true.b)

pp = plot(effs)

pp +
  geom_line(data=r.effs,
            aes(x=x,y=m.eff),linetype=2) + 
  theme_classic()

#####
#calculate overlaps

mm = merge(r.effs,pp$data,by=c('dim','x')) %>%
  mutate(ol=ifelse(m.eff>ul & m.eff<ll,TRUE,FALSE))

