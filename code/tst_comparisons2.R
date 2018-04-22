
###
#calculate R-hat using coda
source('config~.R')

#dependencies
library(apcwin) #custom package in dev
library(lme4) #random effects
library(merTools) #random effects display and sampling
library(parallel) #to speed up sampling
library(dplyr) #data manipulation
library(reshape2) # data manipulation
library(ggplot2) #plotting
library(coda) #mcmc tools

#funciton to generate rhat
rhats = function(sampobj,varname,makeplot=FALSE){
  #var is variable name from summary of sample object
  require(coda)
  

  var = sampobj$summaries[,varname]
  chainlen = length(var)/sampobj$chains
  var = lapply(split(var,ceiling(seq_along(var)/chainlen)),mcmc)
  
  #print(makeplot)
  
  var = as.mcmc.list(var)
  
  if(makeplot){
    print('plotting')
    plot(var)}
  
  return(gelman.diag(var))
  
  }

#rhats(effs$sampobj,'w',makeplot=TRUE)

##################
#estimate implied model; and calculate bayes factor against HAPC

#function to return breaks with a particular sensitivity
estimate_breaks=function(effectsobj,pval=0.05){
  #pval identifies sensitivity for breakpoints

  tdat = effectsobj$sampobj$data
  
  uvals = list(
    a=unique(tdat$a)[order(unique(tdat$a))],
    p=unique(tdat$p)[order(unique(tdat$p))],
    c=unique(tdat$c)[order(unique(tdat$c))]
  )
  
  #fdraw differences from effects object
  dd = tstdiff(effs)
  brk = lapply(dd,function(x){
    res = x$pval<=pval
    #last value is always true for window break
    res = c(res,TRUE)
    return(res)})
  
  
  winbrks = list(
    a=uvals[['a']][brk[['a']]],
    p=uvals[['p']][brk[['p']]],
    c=uvals[['c']][brk[['c']]]
  )
  
  return(winbrks)
  
}


################33
#function to estimate implied models
imply_model = function(effectobj,pval=0.05){
  
  tdat = effectobj$sampobj$data
  
  brks = estimate_breaks(effs,pval=pval)
  brks[['p']] = c(min(tdat$p)-1,brks[['p']])
  brks[['c']] = c(min(tdat$c)-1,brks[['c']])
  
  tdat$ip = window(tdat$p,breaks=brks[['p']])
  tdat$ic = window(tdat$c,breaks=brks[['c']])
  
  return(lm(y~a + I(a^2) + ip + ic,data=tdat))

}

#function to calculate bayes factor
bayesfactor = function(new,old){
  #new and old are models to compare
  
  bic1 = BIC(new); bic2=BIC(old)
  return(exp(bic1-bic2)/2)
  
}

LR = function(new,old){
  
  ll1 = logLik(new); ll2 = logLik(old)
  return(exp(ll1-ll2))
}

#bayesfactor(imod,truemodel)
#bayesfactor(hapc,truemodel)
#bayesfactor(hapc,imod)

#hapc2 = lmer(y~a + I(a^2) + (1|ip) + (1|ic),
#             data=tdat)



############333
#cycle through simulations

#directory

simdir = paste0(outdir,'simdata1/')
resfile = paste0(simdir,'comp.csv')
simtable = read.csv(paste0(simdir,'fullsimtable.csv')) %>%
  arrange(simnumber)

comp = read.csv(paste0(simdir,'comp.csv'))
comp$r.r2 = comp$r.w = 
comp$bf.i001 = comp$bf.i01 = comp$bf.i05 =  
  as.numeric(NA)


sims = simtable$simnumber

for(s in 1:length(sims)){
  print(paste(s,'of',length(sims)))

  lim = comp$simnumber == sims[[s]]  
  #load data
  load(paste0(simdir,'sim',sims[s],'.RData'))
  
  #calculate rhat of r-squared
  comp$r.r2[lim] = rhats(effs$sampobj,'r2')[[1]][1]
  comp$w.r2[lim] = rhats(effs$sampobj,'w')[[1]][1]
  
  #generate implied models at 0.05, 0.01, and 0.001
  i001 = imply_model(effs,pval=0.001)
  i01 = imply_model(effs,pval=0.01)
  i05 = imply_model(effs,pval=0.05)
  
  #is first better than 2d
  comp$bf.i05[lim] = sig(bayesfactor(i05,hapc))
  comp$bf.i01[lim] = sig(bayesfactor(i01,hapc))
  comp$bf.i001[lim] = sig(bayesfactor(i001,hapc))
  
  
}
  
##update comp
ncomp = comp[226:nrow(comp),]  #because I accidentally duplicated
ncomp = merge(ncomp,simtable,by='simnumber')

write.csv(ncomp,paste0(simdir,'ncomp.csv'),row.names=FALSE)

ncomp = ncomp %>% 
  mutate_at(vars(bf.i05,bf.i01,bf.i001),
            funs(
              ifelse(grepl('\\*',.),1,0)
            ))

prop.table(table(ncomp$bf.i001))

#############
#lower & random breaks
summary(glm(bf.i001~p.blocks+c.blocks+p.type+c.type,
    data=ncomp,
    family=binomial('logit')),
    contrasts = list(
      p.type = 'contr.sum',
      c.type = 'contr.sum'
    ))

#summary(glm(bf.i001~factor(p.autocorr)+factor(c.autocorr),
#            data=ncomp,
#            family=binomial('logit')))

############33
#write a table

simsum = cbind(
  ncomp %>%
  group_by(p.blocks,p.type) %>%
  summarize(prop_aws=mean(bf.i001),
            simulations=n()),
  ncomp %>%
    group_by(c.blocks,c.type) %>%
    summarize(prop_aws=mean(bf.i001),
              simulations=n()))

library(knitr)
kable(simsum,digits=2)

