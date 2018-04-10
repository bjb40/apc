

rm(list=ls())
source('config~.R')
library(apcwin)
library(reshape2)
library(msm)

#lnum=34 #best one for hapc
lnum =15 #15 works fine too!

load(paste0(outdir,'simdata_hapc/sim',lnum,'.RData'))
#load(paste0(outdir,'simdata1/sim',lnum,'.RData'))


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
            data = tdat,
            contrasts=list(
              p.hf = 'contr.sum',
              c.hf = 'contr.sum'
            ))

#####################3
#predict

hapc.effs = function(mermod,newdat){
  #mermod is fitted model, newdat is new data
  
  
  
  effs=simulate(mermod,
                nsim=100,
                use.u = TRUE,
                newdat=newdat,
                allow.new.levels=TRUE)

  
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


#set contrasts option; see http://faculty.nps.edu/sebuttre/home/R/contrasts.html

tdat$yc = tdat$y - mean(tdat$y)
tt = lm(yc~ca+ca2+pf+cf,data=tdat,
        contrasts=list(pf='contr.sum',
                       cf = 'contr.sum'))
true.b =list()
beta = coef(tt)

#rename columns of beta ...


###a effects -- note max is for the "omitted" category

xh= data.frame(a=min(tdat$a):max(tdat$a),
               pf = window(mean(tdat$p),breaks=attr(tdat$pf,'breaks')),
               cf = window(mean(tdat$c),breaks=attr(tdat$cf,'breaks'))) 

xh$ca = xh$a-mean(tdat$a); xh$ca2 = xh$ca^2
xm = model.matrix(~ca + ca2 + pf + cf,
                  data=xh,
                  contrasts.arg = list(cf = 'contr.sum',
                                       pf = 'contr.sum'))
facs = grepl('cf|pf',colnames(xm))
xm[,facs] = xm[,facs]*.5


true.b[['a']] = data.frame(
  x=min(tdat$a):max(tdat$a),
  #m.eff = predict(tt,newdata=xh))
  m.eff = xm %*% beta[colnames(xm)])

true.b[['a']]$dim='a'


###p effects
xh= data.frame(ca= 0,
               ca2 = 0,
               pf=scopedummy(tdat$pf,unique.vals=unique(tdat$p)))
xm = model.matrix(~.,data=xh,
                  contrasts.arg=list(pf='contr.sum'))
xm=xm[,2:ncol(xm)]

true.b[['p']] = data.frame(
  x=unique(tdat$p),
  m.eff = xm %*% beta[colnames(xm)])
true.b[['p']]$dim='p'

###c effects

xh= data.frame(ca= 0,
               ca2 = 0,
               cf=scopedummy(tdat$cf))
xm = model.matrix(~.,data=xh,
                  contrasts.arg = list(cf='contr.sum'))
xm = xm[,2:ncol(xm)]
true.b[['c']] = data.frame(
  x=min(tdat$c):max(tdat$c),
  m.eff = xm %*% beta[colnames(xm)])
true.b[['c']]$dim='c'

#####
#plot

r.effs = do.call(rbind,true.b)

#r.effs = r.effs %>% 
#  mutate(m.eff = ifelse(dim=='a',m.eff+21,m.eff))

#redraw effects... using effect coding
effs = draw_sumeffs(effs$sampobj,tol=0.001)
pp = plot(effs)

ages = merge(pp$data %>% filter(dim=='a'),
             r.effs %>% filter(dim=='a'),by='x')

plot(ages[,c('Fit','m.eff')])

##plot comparisons
partp =  pp +
   geom_line(data=r.effs,
              aes(x=x,y=m.eff),linetype=2)

fullp = partp + 
   geom_line(data=pred,
             aes(x=x,y=fit),linetype=3) +
   geom_ribbon(data = pred,
               aes(x=x,y=fit,ymax=ul,ymin=ll),alpha=0.1) +
    theme_classic()

print(fullp)


#real age
rage = simulated$effects[1:2]

#hapc age
hage = coef(summary(hapc))[2:3,]

#estimate of age from simulated posterior


agesim = as.data.frame(t(effs$effects$a))
agesim$a = as.numeric(colnames(effs$effects$a))
agesim = melt(agesim,id='a') %>%
  rename(y=value) %>%
  mutate(ca = a - mean(a),
         ca2 = ca^2) %>%
  dplyr::select(-variable)

wage = coef(summary(lm(y~ca+ca2,data=agesim)))


print(rage)
print(hage)
print(wage)


###########3
#nonmarginal

#nmeffs = draw_effs(effs$sampobj,marginal=FALSE,tol=0.01)

#plot(nmeffs) + geom_line(data=r.effs,
#          aes(x=x,y=m.eff),linetype=2) + 
#  geom_line(data=pred,
#            aes(x=x,y=fit),linetype=3) +
#  geom_ribbon(data = pred,
#              aes(x=x,y=fit,ymax=ul,ymin=ll),alpha=0.1) +
#  theme_classic()

