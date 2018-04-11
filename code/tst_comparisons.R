

rm(list=ls())
source('config~.R')
library(apcwin)
library(reshape2)
library(msm)

#lnum=34 #best one for hapc
lnum =24 #15 works fine too for hapc!

#load(paste0(outdir,'simdata_hapc/sim',lnum,'.RData'))
load(paste0(outdir,'simdata1/sim',lnum,'.RData'))


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

#because of strucure there are some unused year factors --- these are deleted, here
tdat$pf = droplevels(tdat$pf)
#mean center age
tdat$yc = tdat$y - mean(tdat$y)

#recover breaks
class(tdat$pf) = c(class(tdat$pf),'window')
attr(tdat$pf,'breaks') = c(min(tdat$p)-1,
                           simulated$p.breaks[simulated$p.breaks %in% unique(tdat$p)])

hapc = lmer(yc~ca + ca2 +
              p.hf + c.hf +
              (1|p.hf) + (1|c.hf),
            data = tdat,
            contrasts=list(
              p.hf = 'contr.sum',
              c.hf = 'contr.sum'
            ))

convergence = hapc@optinfo$conv$lme4$code

#####################3
#predict -- HAPC

beta = coef(summary(hapc))[,'Estimate']
vc = vcov(hapc)

mna = mean(tdat$ca)
mnp = window(mean(tdat$p),breaks=attr(tdat$p.hf,'breaks'))
mnc = window(mean(tdat$c),breaks=attr(tdat$c.hf,'breaks'))

###a effects -- note max is for the "omitted" category
p.marginals = prop.table(table(tdat$p.hf))
names(p.marginals) = paste0('p.hf',1:length(levels(tdat$p.hf)))
c.marginals = prop.table(table(tdat$c.hf))
names(c.marginals) = paste0('c.hf',1:length(levels(tdat$c.hf)))

marginals=c(p.marginals,c.marginals)

xh= data.frame(a=min(tdat$a):max(tdat$a),
               p.hf = mnp,
               c.hf = mnc) 

marginals = do.call(rbind,lapply(1:nrow(xh),FUN=function(x){return(marginals)}))
marginals = marginals/2

xh$ca = xh$a-mean(tdat$a); xh$ca2 = xh$ca^2
xm = model.matrix(~ca + ca2 + p.hf + c.hf,
                  data=xh,
                  contrasts.arg = list(c.hf = 'contr.sum',
                                       p.hf = 'contr.sum'))
colorder = colnames(xm) #preserve order
facs = grepl('c.hf|p.hf',colnames(xm))
xm= cbind(xm[,!facs],marginals[,colorder[facs]])

pred = list()

pred[['a']] = data.frame(
  x=min(tdat$a):max(tdat$a),
  #m.eff = predict(tt,newdata=xh))
  m.eff = xm %*% beta[colnames(xm)],
  se= sqrt(diag(xm %*% vc[colnames(xm),colnames(xm)] %*% t(xm))))

pred[['a']]$dim = 'a'

###p effects
xh= data.frame(ca= 0,
               ca2 = 0,
               p.hf=scopedummy(tdat$p.hf,unique.vals=unique(tdat$p)))
xm = model.matrix(~.,data=xh,
                  contrasts.arg=list(p.hf='contr.sum'))
xm=xm[,2:ncol(xm)]

pred[['p']] = data.frame(
  x=unique(tdat$p),
  m.eff = xm %*% beta[colnames(xm)],
  se= diag(xm %*% vc[colnames(xm),colnames(xm)] %*% t(xm))
  )
pred[['p']]$dim='p'

###c effects
xh= data.frame(ca= 0,
               ca2 = 0,
               c.hf=scopedummy(tdat$c.hf))
xm = model.matrix(~.,data=xh,
                  contrasts.arg = list(c.hf='contr.sum'))
xm = xm[,2:ncol(xm)]
pred[['c']] = data.frame(
  x=min(tdat$c):max(tdat$c),
  m.eff = xm %*% beta[colnames(xm)],
  se = diag(xm %*% vc[colnames(xm),colnames(xm)] %*% t(xm))
  )
pred[['c']]$dim='c'

pred = do.call(rbind,pred) %>%
  arrange(dim) %>%
  mutate(dim_f = factor(dim, levels=c('a','p','c'),
                        labels=c('Age','Period','Cohort')))

##########3
#real effects

#set contrasts option; see http://faculty.nps.edu/sebuttre/home/R/contrasts.html

tt = lm(yc~ca+ca2+pf+cf,data=tdat,
        contrasts=list(pf='contr.sum',
                       cf = 'contr.sum'))
true.b =list()
beta = coef(tt)

#rename columns of beta ...

###a effects -- note max is for the "omitted" category
p.marginals = prop.table(table(tdat$pf))
names(p.marginals) = paste0('pf',1:length(levels(tdat$pf)))
c.marginals = prop.table(table(tdat$cf))
names(c.marginals) = paste0('cf',1:length(levels(tdat$cf)))

marginals=c(p.marginals,c.marginals)

xh= data.frame(a=min(tdat$a):max(tdat$a),
               pf = window(mean(tdat$p),breaks=attr(tdat$pf,'breaks')),
               cf = window(mean(tdat$c),breaks=attr(tdat$cf,'breaks'))) 

marginals = do.call(rbind,lapply(1:nrow(xh),FUN=function(x){return(marginals)}))
marginals = marginals/2

xh$ca = xh$a-mean(tdat$a); xh$ca2 = xh$ca^2
xm = model.matrix(~ca + ca2 + pf + cf,
                  data=xh,
                  contrasts.arg = list(cf = 'contr.sum',
                                       pf = 'contr.sum'))
colorder = colnames(xm) #preserve order
facs = grepl('cf|pf',colnames(xm))
xm= cbind(xm[,!facs],marginals[,colorder[facs]])


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

r.effs = do.call(rbind,true.b) %>%
  arrange(dim) %>%
  mutate(dim_f = factor(dim, levels=c('a','p','c'),
                        labels=c('Age','Period','Cohort')))

effs = draw_sumeffs(effs$sampobj,tol=0.05)

pp = plot(effs)

##plot comparisons
partp =  pp +
   geom_line(data=r.effs,
              aes(x=x,y=m.eff),size=1.1,color='blue',alpha=0.3)

partp 

hapc_ribbon =    geom_ribbon(data = pred,
                             aes(x=x,y=m.eff,
                                 ymax=m.eff+1.96*se,
                                 ymin=m.eff-1.96*se),alpha=0.1) 

fullp = partp + 
   geom_line(data=pred,
             aes(x=x,y=m.eff),linetype=3) +
    theme_classic()

print(fullp + hapc_ribbon)

#collect differences (i.e. "S", and test against window breaks...)
s.a = do.call(rbind,tstdiff(effs))
attr(tdat$pf,'breaks')
attr(tdat$cf,'breaks')

#true differences in smae dataframe
tdelta = r.effs %>% group_by(dim) %>%
  mutate(real.delta=lead(m.eff) - (m.eff),
         dfname = paste0(dim,'.',lead(x),'-',x)) %>%
  select(real.delta,dfname)

#tests
tdelta = merge(tdelta,s.a,by.x='dfname',by.y='row.names')
quantile(tdelta$real.delta-tdelta$diff,probs=0.5)
quantile(tdelta$real.delta-tdelta$diff,probs=c(0.975,0.025))
sum(tdelta$real.delta-tdelta$diff>0)/nrow(tdelta)

