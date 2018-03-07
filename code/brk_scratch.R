#####
#from analyze gss_parsol

#load recessions data


###########
#deltas

df.cohort = df %>% filter(type=='c')
yrs=df.cohort$id
recessions=read.csv('H:/projects/proposal/r_study/output/recessionlist.csv')
rc = read.csv('H:/projects/proposal/r_study/output/recession.csv')

delts = function(e){
  #calculate proportion >0 across cohorts
  res = e$c[1,]
  res[1] =1
  names(res) = yrs
  
  for(i in 2:ncol(e$c)){
    res[i] = sum(e$c[,i] >e$c[,i-1])/nrow(e$c) #worse off
  }
#  print(e$w)
  res = c(res,e$w)
  names(res) = c(yrs,'wt')
  return(res)
  
}

props = as.data.frame(do.call(rbind,lapply(effects,delts)))
#calculate weighted mean
pp = data.frame(birthyear=yrs,
                prop=apply(props,2,function(x) mean(x)*x[72])[1:71])

pp$propd = abs(pp$prop-.5)

pp = merge(pp %>% filter(birthyear>=1917) 
           ,recessions,by='birthyear',all.x=TRUE,all.y=FALSE)

pp$recession_exp  = factor(pp$recession_exp)

#want to redo for a bayesian p-value in the 
pp %>% group_by(recession_exp) %>% summarize(mean(propd)) 
pp %>% group_by(recession_exp) %>% summarize(mean(prop))

pp$shade = ifelse(pp$prop<=0.025 | pp$prop>=0.975,1,0.6)
tt=table(pp[,c('recession_exp','shade')])
prop.table(tt,margin=1)

ggplot(pp,
       aes(x=birthyear,
           y=prop)) + 
  geom_point(aes(color=recession_exp,shape=recession_exp,alpha=propd*2), 
             size=3) + 
  geom_hline(aes(yintercept=0.025)) +
  geom_hline(aes(yintercept=0.975)) +
  geom_line(alpha=0.5,linetype=2) +
  theme_classic()
  

#pp = lapply(yrs,function(x)
#  wtmn(props,yrs,))

##############
#breaks style -- this is in line with variable selection in BMA

brk = recessions %>% 
  filter(birthyear %in% yrs)
brk$nbreak=0

for(i in seq_along(breaks$c)){
  #print(as.numeric(brk$birthyear %in% breaks$c[i]))
  #calculate weighted proportion
  brk$nbreak = brk$nbreak + 
    as.numeric(brk$birthyear %in% breaks$c[[i]])*win$w_prime[i]
}

brk$prop = brk$nbreak

brk$recession_exp=factor(brk$recession_exp)

ggplot(brk %>% filter(birthyear>=1917),
       aes(x=birthyear,
           y=prop,
           color=recession_exp,
           shape=recession_exp)) + 
  geom_point(size=3)

brk %>% group_by(recession_exp) %>% summarize(mean(prop))

