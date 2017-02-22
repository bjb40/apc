##https://www.r-bloggers.com/dirichlet-process-infinite-mixture-models-and-clustering/

par(mfrow=c(2,2))

#higher alpha leads to higher probability of lower values
stick_draw = function(num.vals, alpha) {
  betas = rbeta(num.vals, 1, alpha)
  stick.to.right = c(1, cumprod(1 - betas))[1:num.vals]
  weights = stick.to.right * betas
  weights
}

#low alpha is low pooling; high alpha is diffuse pooling
#you could invert it with a low alpha
len = 20 #length

winnum=matrix(0,1000,1)
cum.breaks=matrix(0,1000,len)
#first number is windows; 
#second is alpha: sets break points...

for(i in 1:1000){
  d=stick_draw(len,20)
  breaks=unique(c(0,cumsum(round(d*len)),len))
  cum.breaks[i,1:length(breaks)]=breaks
  #  print(breaks)
  a=1:len
  testwin=window(a,breaks=breaks)
  winnum[i] = length(levels(testwin))
}

hist(winnum); table(winnum); mean(winnum)
#print('mean')
#apply(cum.breaks,2,mean,na.rm=TRUE)
#print('max')
#apply(cum.breaks[,1:len-1],2,max,na.rm=TRUE)
#print('min')
#apply(cum.breaks[,1:len-1],2,min,na.rm=TRUE)

library(MCMCpack)

#alpha=rep(1,20)
#symmetric is where alpha is a single scalar
alpha=c(.1,.1,1,1,
        2,2,2,2,
        3,3,3,3,
        10,10,10,10)

d=rdirichlet(1000,alpha)
hist(apply(d,1,function(x) sum(x>(1/20)))) #can use the same function as below
print(colSums(d)/nrow(d))

#unique algorithm...
#need proof???
#

####uniform equality partition ....
#50% probability it is equal to next value
#only about 200 combos
winnum=matrix(0,1000,1)
cum.breaks=matrix(as.numeric(NA),1000,20)

for(i in 1:1000){
  mean.wins=round(runif(1)*20); #print(mean.wins)
  a=1:20
  winprob=1-(mean.wins/max(a))
  partition=runif(19)>winprob
  if(!any(partition)){next} #skip cases with no breaks (could add as continuous...)
  breaks=which(partition==TRUE)
  cum.breaks[i,1:length(breaks)]=breaks
  testwin=window(a,breaks=unique(c(0,breaks,max(a))))
  print(breaks)
  winnum[i] = length(levels(testwin))
  
}

hist(winnum); table(winnum); mean(winnum)
cat('mean window length', 20/mean(winnum))

sum(!duplicated(cum.breaks),margin=1)

print('mean')
apply(cum.breaks,2,mean,na.rm=TRUE)
print('max')
apply(cum.breaks[,1:19],2,max,na.rm=TRUE)
print('min')
apply(cum.breaks[,1:19],2,min,na.rm=TRUE)




####alt idea ---- each is an ind parameter of equality to next value
#then in the mh algorithm, they can be updated
#you can probably figure out using combinotorics the best way 
# to limit window size
# e.g. if probability of equality is .2 ... sum(.2>runif(1000))/1000

#reasonable prior: as likely to equal as not equal...
#should look at that bayesian model from the green book some time...

winnum=matrix(0,1000,1)
cum.breaks=matrix(as.numeric(NA),1000,20)


p=rep(.5,20) # can draw from a proposal support between 0 and 1
#can also use this to sample breaks --  looks really close do dirichelet

for(i in 1:1000){
  equ.draws = runif(20)
  partition = p>equ.draws
  if(!any(partition)){next} #skip cases with no breaks (could add as continuous...)
  breaks=c(0,which(partition==TRUE),20)
  print(breaks)  
  cum.breaks[i,1:length(breaks)]=breaks
  testwin=window(a,breaks=unique(breaks))
  print(breaks)
  winnum[i] = length(levels(testwin))
  
}

hist(winnum); table(winnum); mean(winnum)
cat('mean window length', 20/mean(winnum))


print('mean')
apply(cum.breaks,2,mean,na.rm=TRUE)
print('max')
apply(cum.breaks[,1:19],2,max,na.rm=TRUE)
print('min')
apply(cum.breaks[,1:19],2,min,na.rm=TRUE)



####
#equality constraints sampler; i.e. sampling "breaks" directly
#doesn't work
winnum=matrix(0,1000,1)
cum.breaks=matrix(as.numeric(NA),1000,20)

n.constraints=10
for(i in 1:1000){
  #set uniform partition
  partition = rep(TRUE,length(a)-1)
  #draw random constraints based on n.constraints, 
  #i.e. 1 constraint = equal to the next beta (behind may be more natural)
  partition[sample(head(a,-1),n.constraints)] = FALSE
  if(!any(partition)){next} #skip cases with no breaks (could add as continuous...)
  breaks=c(0,which(partition==TRUE),20)
  print(breaks)  
  cum.breaks[i,1:length(breaks)]=breaks
  testwin=window(a,breaks=unique(breaks))
  print(breaks)
  winnum[i] = length(levels(testwin))
  
}

hist(winnum); table(winnum); mean(winnum)
cat('mean window length', 20/mean(winnum))


print('mean')
apply(cum.breaks,2,mean,na.rm=TRUE)
print('max')
apply(cum.breaks[,1:19],2,max,na.rm=TRUE)
print('min')
apply(cum.breaks[,1:19],2,min,na.rm=TRUE)

sum(!duplicated(cum.breaks),margin=1)
