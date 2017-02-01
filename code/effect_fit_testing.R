##42 is single windows for a and p with no c y4 
##7 is single windows for p and c with no a y2

testmodel=best
actual='s1'
dv='y1'
totest = c('a','p','c')

#one off for now based on looking...
xhat=data.frame(a=window(tdat$a,win[best,'a']),
                p=window(tdat$p,win[best,'p']),
                c=window(tdat$c,win[best,'c']))
xhat=model.matrix(~.,xhat)

#apc_lm(y4~a+p,tdat,age='a',per='p')

for(d in totest){
  
  preds[[d]]=as.data.frame(apply(effects[[testmodel]][[d]],2,mean))
  colnames(preds[[d]])='est'
  rng=apply(effects[[testmodel]][[d]],2,quantile,c(0.025,0.975))
  preds[[d]]$up = rng[2,]
  preds[[d]]$down = rng[1,]
  #s1-s4 are for scenarios --- needs to match with y1-y4
  #preds[[d]]$actual=pltdat[[d]]$s1[order(pltdat[[d]]$id)]
  preds[[d]]$actual=pltdat[[d]][,actual][order(pltdat[[d]]$id)]
  preds[[d]]$id=pltdat[[d]]$id[order(pltdat[[d]]$id)]
  
  #these are the margins / se is pretty close
  
  #need to change axes and limits...actually better to melt into df 
  #and not to use grid.arrange... (same for below)
  calcplt=ggplot(preds[[d]],
                 aes_string(y='est',x='id') ) + 
    geom_point() + 
    geom_errorbar(aes(ymax=up,ymin=down)) +
    geom_line(aes(y=actual))
  
  #print(calcplt)
  #Sys.sleep(2)
  
  best.plt[[d]]=calcplt
  
}


grid.arrange(best.plt[['a']],
             best.plt[['p']],
             best.plt[['c']],
             ncol=2)

#Sys.sleep(5)


#@@@@@@posterior predictive and bayesian p-values

#simulate y distributions from estimates

ytilde = xhat %*% t(as.matrix(allmods[[testmodel]]$betas))

ytilde = sapply(1:1000,function(i)
              ytilde[,i] + rnorm(10000,mean=0,sd=allmods[[testmodel]]$sigma[i]))

#omnibus bayesian pvalues
print(summary(tdat$y1))
print(summary(as.vector(ytilde)))

#mean
sum(apply(ytilde,2,mean)<mean(tdat$y1))/ncol(ytilde)
#sum(apply(ytilde,2,max)<max(tdat$y2))/ncol(ytilde)
#sum(apply(ytilde,2,min)>min(tdat$y2))/ncol(ytilde)


#bayesian p-values by period 
library(dplyr)

actual.period = aggregate(tdat[,dv],by=list(tdat$p),mean)
ytilde.period =  sapply(1:ncol(ytilde),function(i)
                        aggregate(ytilde[,i],by=list(tdat$p),mean)[[2]]
                        )
sapply(seq_along(actual.period[[1]]),function(i)
              sum(ytilde.period[i,]>actual.period[i,2])/ncol(ytilde)
)


#bayesian p-values by cohort 
actual.cohort = aggregate(tdat[,dv],by=list(tdat$c),mean)
ytilde.cohort =  sapply(1:ncol(ytilde),function(i)
  aggregate(ytilde[,i],by=list(tdat$c),mean)[[2]]
)
sapply(seq_along(actual.cohort[[1]]),function(i)
  sum(ytilde.cohort[i,]>actual.cohort[i,2])/ncol(ytilde)
)


#bayesian p-values by age 
actual.age = aggregate(tdat[,dv],by=list(tdat$a),mean)
ytilde.age =  sapply(1:ncol(ytilde),function(i)
  aggregate(ytilde[,i],by=list(tdat$a),mean)[[2]]
)
sapply(seq_along(actual.age[[1]]),function(i)
  sum(ytilde.age[i,]>actual.age[i,2])/ncol(ytilde)
)

par(mfrow=c(1,3))
plot(actual.age); lines(actual.age[[1]],apply(ytilde.age,1,mean))
plot(actual.period); lines(actual.period[[1]],apply(ytilde.period,1,mean))
plot(actual.cohort); lines(actual.cohort[[1]],apply(ytilde.cohort,1,mean))
