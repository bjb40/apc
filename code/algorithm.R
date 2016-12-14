






####from before



#plot facets for drawn data of four scenarios (loess smoothing)

#a = reshape(dat,idvar=)

d=ggplot(dats[[1]],aes(x=a)) + 
  geom_smooth(aes(y=y1,col='senario 1')) +
  geom_smooth(aes(y=y2,col='senario 2')) +
  geom_smooth(aes(y=y3,col='senario 3')) +
  geom_smooth(aes(y=y4,col='senario 4')) 

print(d)

#id time for OLS estimation of single window constraint
t=Sys.time()

#cut window constraints (per luo and hodges)
x = dats[[1]] %>%
  mutate(
    a=factor(a),
    p=cut_interval(p,length=2),
    c=cut_interval(c,length=5)
  ) %>%
  select(a,p,c,y1,y2,y3,y4)

e.beta = list(
  lm(y1~a+p+c,data=x),
  lm(y2~a+p+c,data=x),
  lm(y3~a+p+c,data=x),
  lm(y4~a+p+c,data=x)
  
)

delta.t = Sys.time()-t

cat('One OLS analysis of a single instance window constraint took:', delta.t,'seconds \n\n')
e.time=100*delta.t
cat('Estimated time for one set of constraints:',e.time/60,'minutes \n\n')

#even bins only (with remainder) --- ignoring double cohort//need function... 
configs=10^3 #need to really think about this... windows are adjacent.../if I make unequal
cat('Estimated time for',configs,'window constraints',(e.time*configs)/60/60,'hours\n\n')


#some extra crap re ... compoosition by age (Raudenbush-style)
x$age = as.numeric(levels(x$a)[x$a])

print(
  x %>% group_by(p) %>% summarize(mean(age))
)

print(
  x %>% group_by(c) %>% summarize(mean(age))
)

x1 = x %>% group_by(c) %>% mutate(., mcage = age - mean(age),cage=mean(age))

x1 %>% group_by(c) %>% summarize(reg=mean(age),mc=mean(mcage))

print(
  summary(lm(y1~a+p+c-1,dat=x1))
)

print(
  summary(lm(y1~p+c*mcage+mcage^2-1,dat=x1))
)

####ind v. corr

est_windows = function(){
  
}

comp.dats = list(uncorrelated=inddats[[1]],correlated=cordats[[1]])
for(d in 1:2){
  colnames(comp.dats[[d]])[5:6]= c('c','c2')
}

c.beta = lapply(comp.dats,FUN=function(x) 
  lm(y1~a+a2+p+p2+c+c2,dat=x,x=TRUE))