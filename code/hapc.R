#############
#apc models

rm(list=ls())
library(apcwin)
library(lme4)

source('config~.R')

simdir = paste0(outdir,'simdata_hapc/')
#simdir = paste0(outdir,'simdata1/')


#collect overview of models
simtable = read.csv(
  paste0(simdir,'simtable.csv')
  )

for(lnum in simtable$simnumber){
#for(lnum in 4:4){

rownum = which(simtable$simnumber == lnum)

####load
load(paste0(simdir,'sim',lnum,'.RData'))

tdat = simulated$simdat
tdat$y = as.numeric(tdat$y)

library(lme4)
library(dplyr)

#center age

tdat = tdat %>% 
  mutate(ca = a-mean(a,na.rm=TRUE),
         ca2 = ca^2,
         #real breaks
         pf = window(p,breaks=simulated$p.breaks),
         cf = window(c,breaks=simulated$c.breaks),
         #hapc breaks
         p.hf = window(p,winlength=5),
         c.hf = window(c,winlength=3))

hapc = lmer(y~ca + I(ca^2) +
              p.hf + c.hf +
              (1|p.hf) + (1|c.hf),
            data = tdat)

#####################3
#predict

hapc.effs = function(mermod,newdat){
  #mermod is fitted model, newdat is new data
  
  effs=simulate(mermod,
                nsim=100,
                use.u = TRUE,
                newdat=newdat,
                allow.new.levels=TRUE,
                se.fit=TRUE)
  
  pred = data.frame(fit=apply(effs,1,mean),
                    ll = apply(effs,1,quantile,prob=0.025),
                    ul = apply(effs,1,quantile,prob=0.975))
  
  return(pred)
}


mna = mean(tdat$ca)
mnp = window(mean(tdat$p),breaks=attr(tdat$p.hf,'breaks'))
mnc = window(mean(tdat$c),breaks=attr(tdat$c.hf,'breaks'))

xm.a = data.frame(ca = min(tdat$ca):max(tdat$ca),
                  p.hf=mnp,
                  c.hf=mnc
                  )

pred = hapc.effs(hapc,xm.a) %>%
  mutate(x = unique(tdat$a),
         dim= 'a')

xm.p = data.frame(ca = mna,
                  p.hf=scopedummy(tdat$p.hf,unique.vals=unique(tdat$p)),
                  c.hf=mnc)

pred = rbind(pred,hapc.effs(hapc,xm.p) %>% 
               mutate(x=unique(tdat$p),
                      dim='p'))


xm.c = data.frame(ca = mna,
                  p.hf=mnp,
                  c.hf=scopedummy(tdat$c.hf))

pred = rbind(pred,hapc.effs(hapc,xm.c) %>%
               mutate(x=unique(tdat$c),
                      dim='c'))

##########3
#real effects


#model is identifiable!! (can check directly from betas..., but using scopedummy and predict)
tdat$a2 = tdat$a^2
tt = lm(y~a+a2+pf+cf,data=tdat)
true.b =list()

###a effects
xh= data.frame(a=min(tdat$a):max(tdat$a)); xh$a2 = xh$a^2
xh$pf = window(mean(tdat$p),breaks=attr(tdat$pf,'breaks'))
xh$cf = window(mean(tdat$c),breaks=attr(tdat$cf,'breaks'))

true.b[['a']] = data.frame(
  x=xh$a,
  m.eff = predict(tt,newdata=xh))
true.b[['a']]$dim='a'

###p effects
xh= data.frame(pf=window(unique(tdat$p),
                         breaks=attr(tdat$pf,'breaks')))
xh$a = mean(tdat$a); xh$a2=mean(tdat$a2)
xh$cf = window(mean(tdat$c),breaks=attr(tdat$cf,'breaks'))

true.b[['p']] = data.frame(
  x=unique(tdat$p),
  m.eff = predict(tt,newdata=xh))
true.b[['p']]$dim='p'

###c effects
xh= data.frame(cf=scopedummy(tdat$cf))
xh$a = mean(tdat$a); xh$a2=mean(tdat$a2)
xh$pf = window(mean(tdat$p),breaks=attr(tdat$pf,'breaks'))

true.b[['c']] = data.frame(
  x=min(tdat$c):max(tdat$c),
  m.eff = predict(tt,newdata=xh))
true.b[['c']]$dim='c'

#####
#plot

r.effs = do.call(rbind,true.b)


############
#plot (and calcualate overlaps)



pp = plot(effs)

#plot comparisons
#pp +
#   geom_line(data=r.effs,
#             aes(x=x,y=m.eff),linetype=2) + 
#   geom_line(data=pred,
#             aes(x=x,y=fit),linetype=3) +
#   geom_ribbon(data = pred,
#               aes(x=x,y=fit,ymax=ul,ymin=ll),alpha=0.1) +
#    theme_classic()

##########3
#calculate differences
mm = merge(r.effs,pp$data,by=c('dim','x')) %>%
  mutate(ol=ifelse(m.eff>ul & m.eff<ll,TRUE,FALSE))

mm = mm %>%
  rename(ll = ul,
         ul = ll)

pred = pred %>% 
  rename(Fit.h = fit,
         ll.h = ll,
         ul.h = ul)

mm = merge(mm,pred,by=c('dim','x')) %>%
  mutate(ol.h = ifelse(m.eff>ll.h & m.eff<ul.h,TRUE,FALSE),
         rng.h = abs(ul.h-ll.h),
         rng = abs(ul-ll))

#hapc results
simtable[rownum,'hol'] = mean(mm$ol.h)
simtable[rownum,'hrmse'] = summary(hapc)$sigma
simtable[rownum,'efficiency'] = mean(mm$rng/mm$rng.h)

#save simtable results
write.csv(simtable,
          file=paste0(simdir,'fullsimtable.csv'),
          row.names=FALSE)
}

