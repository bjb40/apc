##https://www.r-bloggers.com/dirichlet-process-infinite-mixture-models-and-clustering/


#higher alpha leads to higher probability of lower values
stick_draw = function(num.vals, alpha) {
  betas = rbeta(num.vals, 1, alpha)
  stick.to.right = c(1, cumprod(1 - betas))[1:num.vals]
  weights = stick.to.right * betas
  weights
}

winnum=matrix(0,1000,1)
cum.breaks=matrix(0,1000,20)
#first number is windows; 
#second is alpha: sets break points...

for(i in 1:1000){
  d=stick_draw(20,5)
  breaks=unique(c(0,cumsum(round(d*20)),20))
  cum.breaks[i,1:length(breaks)]=breaks
  #  print(breaks)
  a=1:20
  testwin=window(a,breaks=breaks)
  winnum[i] = length(levels(testwin))
}

hist(winnum); table(winnum); mean(winnum)
print('mean')
apply(cum.breaks,2,mean,na.rm=TRUE)
print('max')
apply(cum.breaks[,1:19],2,max,na.rm=TRUE)
print('min')
apply(cum.breaks[,1:19],2,min,na.rm=TRUE)


#unique algorithm...
#need proof???
#

####uniform equality partition ....
#50% probability it is equal to next value
winnum=matrix(0,1000,1)


cum.breaks=matrix(as.numeric(NA),1000,20)

for(i in 1:1000){
  mean.wins=round(runif(1)*20); #print(mean.wins)
  a=1:20
  winprob=1-(mean.wins/max(a))
  partition=runif(19)>winprob
  if(!any(partition)){next} #skipp cases with no breaks (could add as continuous...)
  breaks=which(partition==TRUE)
  cum.breaks[i,1:length(breaks)]=breaks
  testwin=window(a,breaks=c(0,breaks,max(a)))
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
