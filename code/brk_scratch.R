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
  

#add 18 and 22 to get a sense of 
pp = pp %>%
  mutate(year18 = birthyear+18,
         year22 = birthyear+22)

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


#############
#pptable

rc$minbirth =rc$yr_begin-22
rc$maxbirth =rc$yr_end-18

View(rc %>% filter(yr_end>1908) %>% arrange(length))
View(pp %>% select(birthyear,prop,year18,year22) %>% arrange(prop) %>% mutate(round(prop,2)))
##note wwii = Dec 1941 (1942) to 1945

#ggplot(rc %>% filter(minbirth>1917)) +
#  geom_segment(aes(x=minbirth,y=length/12,xend=maxbirth,yend=length/12)) +
#  geom_point(data=pp %>% filter(prop<.05),aes(x=birthyear),y=1,size=4)

####this one, but you want the size of the effect....

ggplot(rc %>% filter(minbirth>1917)) +
  geom_segment(aes(x=yr_begin,y=1,xend=yr_end,yend=1),size=4) +
  geom_point(data=pp %>% filter(prop<.05),aes(x=birthyear+18),y=1,size=4,shape=1) 

ggplot(rc %>% filter(minbirth>1917)) +
  geom_segment(aes(x=yr_begin,y=1,xend=yr_end,yend=1),size=2) +
  geom_point(data=pp %>% filter(prop<.05),aes(x=year22),y=1,size=2,shape=1) 


ggplot(rc %>% filter(minbirth>1917)) +
  geom_segment(data=pp %>% filter(prop<.05),aes(x=year18,y=1,xend=year22,yend=1)) +
  geom_segment(aes(x=yr_begin,y=1,xend=yr_end,yend=1),size=2) 


#ggplot(melt(rc %>% filter(minbirth>1917) %>% select(yr_end,yr_begin,length),id='length')) +
#  geom_point(aes(x=value,y=1,color=variable,size=length/12)) +
#  geom_segment(data=pp %>% filter(prop<.05),aes(x=year18,y=1,xend=year22,yend=1)) 

pp$rlength = 0

for(r in 1:nrow(rc)){
#  byrs = (rc$yr_begin[r]-22):(rc$yr_end[r]-18)
  print(c(r,byrs))
  pp$rlength = ifelse(pp$birthyear %in% byrs,rc$length[r]/12,pp$rlength)
}

pplim = pp %>% filter(prop<.05)

ggplot(pp,
       aes(x=birthyear,
           y=prop)) + 
  geom_point(aes(color=rlength,shape=recession_exp)) + 
  geom_hline(aes(yintercept=0.025)) +
  geom_hline(aes(yintercept=0.975)) +
  geom_line(alpha=0.5,linetype=2) +
  theme_classic()

ggplot(pp,aes(y=prop,x=rlength,alpha=propd*2)) +
         geom_point(size=2) + theme_classic()

ggplot(pp,aes(color=prop<0.05,x=rlength)) +
  geom_bar() + theme_classic()

