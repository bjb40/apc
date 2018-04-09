########
#for new dgp
#note this relies on variables generated in from simulator!


est = apcsamp(dat=tdat,dv='y',
              cores=4,method='ml',
              chains=4,samples=250)

#summary(est)

use.marginal=TRUE

effs = draw_effs(est,tol=0.01,marginal=use.marginal)


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
xh= data.frame(pf=window(unique(tdat$p),
                         breaks=attr(tdat$pf,'breaks')))
xh$a = mean(tdat$a); xh$a2=mean(tdat$a2)
xh$cf = window(mean(tdat$c),breaks=cdraw$breaks)

true.b[['p']] = data.frame(
  x=unique(tdat$p),
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

#pp +
#  geom_line(data=r.effs,
#            aes(x=x,y=m.eff),linetype=2) + 
#  theme_classic()

#####
#calculate overlaps

mm = merge(r.effs,pp$data,by=c('dim','x')) %>%
  mutate(ol=ifelse(m.eff>ul & m.eff<ll,TRUE,FALSE))

##############
#calculate breaks

allb = lapply(effs$sampobj$breaks,function(x)
  x[effs$sampled])

rng=lapply(c('a','p','c'),
           function(x) range(effs$sampobj$data[,x]))
names(rng)=c('a','p','c')

brks = list()
for(d in names(allb)){
  brks[[d]] = lapply(allb[[d]], function(x)
    rng[[d]][1]:rng[[d]][2] %in% x)
  brks[[d]] = do.call(rbind,brks[[d]])
}

ps = colMeans(brks[['p']])
pp = min(tdat$p):max(tdat$p) %in% simulated$p.breaks
p.miss = mean(ps-pp)

cs = colMeans(brks[['c']])
cp = min(tdat$c):max(tdat$c) %in% simulated$c.breaks
c.miss = mean(cs-cp)

exact.miss = 1-mean(cs==cp) #what is there now -- OR
av.miss = mean(abs(cs-cp)) #this is the total "deviance" proportion

