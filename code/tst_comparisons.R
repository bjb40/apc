
#rm(list=ls())
source('config~.R')

#dependencies
library(apcwin) #custom package in dev
library(lme4) #random effects
library(merTools) #random effects display and sampling
library(parallel) #to speed up sampling
library(dplyr) #data manipulation
library(reshape2) # data manipulation
library(ggplot2) #plotting

#lnum =210 #15 works fine too for hapc .. 36 c

#load(paste0(outdir,'simdata_hapc/sim',lnum,'.RData'))
#load(paste0(outdir,'simdata1/sim',lnum,'.RData'))


tdat = simulated$simdat
tdat$y = as.numeric(tdat$y)


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
attr(tdat$pf,'breaks') = simulated$p.breaks

#can rescale; doesn't seem to fix convergence
tdat$cal = tdat$ca/10
tdat$cal2 = tdat$cal^2

#model with factors in r/e translating 
#yang ccrem models from sas:
#http://yangclaireyang.web.unc.edu/age-period-cohort-analysis-new-models-methods-and-empirical-applications/chapter-7/
hapc = lmer(yc~ca + ca2 +
              (1|p.hf) + (1|c.hf),
            data = tdat)

convergence = hapc@optinfo$conv$lme4$code

#####################3
#predict -- HAPC

#already relative to grand mean, 
#because grand mean centered
hapc.re = REsim(hapc)

pred=list()
pred[['p']] = data.frame(
  x=unique(tdat$p)[order(unique(tdat$p))],
  p.hf=scopedummy(tdat$p.hf,unique.vals=unique(tdat$p))
)
pred[['p']] = merge(pred[['p']],hapc.re,by.x='p.hf',by.y='groupID') %>%
  mutate(est = as.numeric(mean),
         ul = est + 1.96*as.numeric(sd),
         ll = mean - 1.96*as.numeric(sd),
         se = as.numeric(sd), #bayesian estimate of se
         dim='p') %>%
  dplyr::select(est,x,dim,ul,ll,se)

#cohort effects
pred[['c']] = data.frame(
  x=unique(tdat$c)[order(unique(tdat$c))],
  c.hf=scopedummy(tdat$c.hf,unique.vals=unique(tdat$c))
)
pred[['c']] = merge(pred[['c']],hapc.re,by.x='c.hf',by.y='groupID') %>%
  mutate(est = as.numeric(mean),
         ul = est + 1.96*as.numeric(sd),
         ll = mean - 1.96*as.numeric(sd),
         se = as.numeric(sd),
         dim='c') %>%
  dplyr::select(est,x,dim,ul,ll,se)

###a effects -- fixed effects -- not held at marginals
#because that is hella hard --- probably simulate the marginals
#from the distributions above (REsim), and generate a 
#distribution .... or can just adjust the "center" based
#on the mean values, and note that they are too conservative
#opting for 2 because easier computationally and theoretically
beta = coef(summary(hapc))[,'Estimate']
vc = vcov(hapc)

p.marginals = prop.table(table(tdat$p.hf))
names(p.marginals) = paste0('p.hf',names(p.marginals))
c.marginals = prop.table(table(tdat$c.hf))
names(c.marginals) = paste0('c.hf',names(c.marginals))

marginals=c(p.marginals,c.marginals)

xh= data.frame(a=min(tdat$a):max(tdat$a)) 

marginals = do.call(rbind,lapply(1:nrow(xh),FUN=function(x){return(marginals)}))
#marginals = marginals/2

xh$ca = xh$a-mean(tdat$a); xh$ca2 = xh$ca^2
xm = model.matrix(~ca + ca2,
                  data=xh)

#standard error -- doesn't account for 2d level uncertainty
#too consevative...
xh$se = sqrt(diag(xm %*% as.matrix(vc) %*% t(xm)))

#estimate adjusted by marginals
beta = c(beta,hapc.re$mean)
names(beta) = c(colnames(vc),
                paste0(hapc.re$groupFctr,hapc.re$groupID))

xm = cbind(xm,marginals)
xh$est = as.vector(xm %*% beta[colnames(xm)])

pred[['a']] = xh %>%
  mutate(x=a,
         ul = est + 1.96*se,
         ll = est - 1.96*se,
         dim='a') %>% 
  dplyr::select(est,x,ul,ll,dim,se)
  


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
vc = vcov(tt)

#rename columns of beta ...

###a effects -- note max is for the "omitted" category
p.marginals = prop.table(table(tdat$pf))
names(p.marginals) = paste0('pf',1:length(levels(tdat$pf)))
c.marginals = prop.table(table(tdat$cf))
names(c.marginals) = paste0('cf',1:length(levels(tdat$cf)))

#contr.sum deletes last category by default
marginals=c(p.marginals[1:length(levels(tdat$pf))-1],
            c.marginals[1:length(levels(tdat$cf))-1])

xh= data.frame(a=min(tdat$a):max(tdat$a),
               pf = window(mean(tdat$p),breaks=attr(tdat$pf,'breaks')),
               cf = window(mean(tdat$c),breaks=attr(tdat$cf,'breaks'))) 

marginals = do.call(rbind,lapply(1:nrow(xh),FUN=function(x){return(marginals)}))
#marginals = marginals/2

xh$ca = xh$a-mean(tdat$a); xh$ca2 = xh$ca^2
xm = model.matrix(~ca + ca2,
                  data=xh)

xm= cbind(xm,marginals)

lm=colnames(xm)
true.b[['a']] = data.frame(
  x=min(tdat$a):max(tdat$a),
  #m.eff = predict(tt,newdata=xh))
  m.eff = xm %*% beta[colnames(xm)],
  se = sqrt(diag(xm %*% vc[lm,lm] %*% t(xm))))

true.b[['a']]$dim='a'


###p effects
xh= data.frame(ca= 0,
               ca2 = 0,
               pf=droplevels(scopedummy(tdat$pf,
                             unique.vals=unique(tdat$p))))

xm = model.matrix(~.,data=xh,
                  contrasts.arg=list(pf='contr.sum'))
xm=xm[,2:ncol(xm)]

lm = colnames(xm)
true.b[['p']] = data.frame(
  x=unique(tdat$p),
  m.eff = xm %*% beta[colnames(xm)],
  se = sqrt(diag(xm %*% vc[lm,lm] %*% t(xm))))
true.b[['p']]$dim='p'

###c effects

xh= data.frame(ca= 0,
               ca2 = 0,
               cf=scopedummy(tdat$cf))
xm = model.matrix(~.,data=xh,
                  contrasts.arg = list(cf='contr.sum'))
xm = xm[,2:ncol(xm)]
lm = colnames(xm)
true.b[['c']] = data.frame(
  x=min(tdat$c):max(tdat$c),
  m.eff = xm %*% beta[colnames(xm)],
  se=sqrt(diag(xm %*% vc[lm,lm] %*% t(xm))))
true.b[['c']]$dim='c'

#####
#plot

r.effs = do.call(rbind,true.b) %>%
  arrange(dim) %>%
  mutate(dim_f = factor(dim, levels=c('a','p','c'),
                        labels=c('Age','Period','Cohort')))

effs = draw_sumeffs(effs$sampobj,
                    tol=0.01)

pp = plot(effs)

#dt=pp$data

"
#commented out while I run the gather
##plot comparisons
th = theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1))

reff_plt = geom_line(data=r.effs,
              aes(x=x,y=m.eff),
             size=1.1,color='blue',alpha=0.3) 

reff_ribbon = 
  geom_ribbon(data=r.effs,
              aes(x=x,y=NULL,
                  ymax=m.eff + 1.96*se,
                  ymin=m.eff-1.96*se),
              fill='blue',alpha=0.15)

hapc_ribbon =   geom_errorbar(data = pred,
                             aes(x=x,y=est,
                                 ymax=ul,
                                 ymin=ll),
                             linetype=3) 

hapc_plt =  geom_line(data=pred,
                 aes(x=x,y=est),linetype=3) 


linecomp = plot(effs,adj.se=FALSE) + hapc_plt + reff_plt
hapc.comp = ggplot() + hapc_plt + hapc_ribbon + reff_plt + reff_ribbon +
  th + facet_wrap(~dim_f,scales='free_x')
sws.comp = pp + reff_plt + reff_ribbon +
  th + facet_wrap(~dim_f,scales='free_x')
cross.comp = pp + hapc_plt + hapc_ribbon
"

#################3
#@test values to save
#################

######
#generate omnibus dataframe with sumamry stats on preds
#use crit??
summdat = pp$data %>%
  rename(sws.est = Fit,
         sws.se = se
         ) %>%
  mutate(sws.ul = sws.est + 1.96*sws.se,
         sws.ll = sws.est - 1.96*sws.se) %>%
  dplyr::select(sws.est,sws.ul,sws.ll,dim,x,dim_f,sws.se)
  
summdat = merge(summdat,pred %>% 
  rename(hapc.est = est,
         hapc.ul = ul,
         hapc.ll = ll,
         hapc.se = se) %>%
  dplyr::select(hapc.est,hapc.ul,hapc.ll,hapc.se,dim,x,dim_f))

summdat =merge(summdat, r.effs %>%
    rename(est = m.eff) %>%
    mutate(ul = est + 1.96*se,
           ll = est - 1.96*se) %>%
  dplyr::select(est,ul,ll,dim,x,dim_f,se))

###
#overlap
#one is for estimate; 2 is for bars
oltest = function(ul1,ll1,ul2,ll2){
  #this tests four conditions for overlap of ci
  #returns logical if there is an overlap
  #overlap is whether 1 overlaps with 2
  #ul1 overlaps
  cond1 = ul1<ul2 & ul1>ll2
  #ll1 overlaps
  cond2 = ll1<ul2 & ll1>ll2
  #nested within
  cond3 = ul1>ul2 & ll1<ll2
  #nested wihtout
  cond4 = ul1<ul2 & ll1>ll2
  eval = do.call(cbind,list(cond1,cond2,cond3,cond4))
  return(as.vector(apply(eval,1,any)))  
  #return(any(c(cond1,cond2,cond3,cond4)))
  
}

#unit tests
#c1 = c(-10,10,5,15,TRUE)
#c2 = c(-10,10,-15,5,TRUE)
#c3 = c(-10,10,-5,5,TRUE)
#c4 = c(-5,5,-10,10,TRUE)
#c5 = c(-20,-10,10,20,FALSE)
#oltests = data.frame(do.call(rbind,list(c1,c2,c3,c4,c5)))
#colnames(oltests) = c('ll1','ul1','ll2','ul2','correct')
#tst = oltests %>%
#  mutate(tested=
#           as.numeric(oltest(ul1,ll1,ul2,ll2)))
#all(tst$tested == tst$correct)

overlap = summdat %>%
  mutate(sws.ol1 = est<sws.ul & est>sws.ll,
         hapc.ol1 = est<hapc.ul & est>hapc.ll,
         sws.ol2 = oltest(ul,ll,sws.ul,sws.ll),
         hapc.ol2 = oltest(ul,ll,hapc.ul,hapc.ll)) %>%
  dplyr::select(sws.ol1,hapc.ol1,sws.ol2,hapc.ol2)


####
#Difference of point estimates (for sign and magnitude)
sm = summdat %>%
  mutate(sws.d = (est-sws.est)/sd(tdat$y),
   hapc.d = (est-hapc.est)/sd(tdat$y),
   sws.s = sign(est) == sign(sws.est),
   hapc.s = sign(est) == sign(hapc.est)
  ) %>% dplyr::select(sws.s,hapc.s,sws.d,hapc.d)

overview = c(colMeans(overlap),colMeans(sm))

"
not sure if these add anything
#collect differences (i.e. 'S', and test against window breaks...)
#i.e. 1 step slope...
s.a = do.call(rbind,tstdiff(effs))
attr(tdat$pf,'breaks')
attr(tdat$cf,'breaks')

#true differences in smae dataframe
tdelta = r.effs %>% group_by(dim) %>%
  mutate(real.delta=lead(m.eff) - (m.eff),
         dfname = paste0(dim,'.',lead(x),'-',x)) %>%
  dplyr::select(real.delta,dfname)

#tests
tdelta = merge(tdelta,s.a,by.x='dfname',by.y='row.names')
quantile(tdelta$real.delta-tdelta$diff,probs=0.5)
quantile(tdelta$real.delta-tdelta$diff,probs=c(0.975,0.025))
sum(tdelta$real.delta-tdelta$diff>0)/nrow(tdelta)
"
#########3
#BIC analysis

sws.bic = mean(effs$fit$bic)
hapc.bic = BIC(hapc)
true.bic = BIC(tt)

BF.sws = exp(sws.bic - true.bic)/2
BF.hapc = exp(hapc.bic - true.bic)/2 
BF = exp(sws.bic - hapc.bic)/2

print('done')
