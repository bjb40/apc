##42 is single windows for a and p with no c y4 
##7 is single windows for p and c with no a y2

testmodel=7
actual='s2'
dv='y2'
totest = c('p','c')

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


grid.arrange(#best.plt[['a']],
             best.plt[['p']],
             best.plt[['c']],
             ncol=2)

Sys.sleep(5)

#x estimates
x = tdat[,totest]
x$p=window(x$p,winlength=1)
x$c=window(x$c,winlength=1)

xmat=model.matrix(~.,data=as.data.frame(x))
#xhat = xhats[[testmodel]]
#u.tdat=unique(tdat[,c('a','p','c')])
u.tdat = tdat[,c('a','p','c')]
xhat= data.frame(p=window(u.tdat$p,1),
                     c=window(u.tdat$c,1))
xhat=model.matrix(~.,xhat)

 summary(lm(tdat[,dv]~xmat-1))

 rev=apply(xhat %*% t(as.matrix(allmods[[testmodel]]$betas)),1,mean)

 #need to weight the mean by its frequency?
 u.tdat$pred=rev
 pred=list(p=aggregate(rev,by=list(u.tdat$p),mean)[,2],
           c=aggregate(rev,by=list(u.tdat$c),mean)[,2])
 
   
#rev=list(p=apply(round(xhat[['p']]) %*% t(as.matrix(allmods[[testmodel]]$betas)),1,mean),
#         c=apply(round(xhat[['c']]) %*% t(as.matrix(allmods[[testmodel]]$betas)),1,mean))

par(mfrow=c(1,2))
plot(1:20,pred[['p']]); lines(1:20,preds[['p']]$actual)
plot(1:39,pred[['c']]); lines(1:39,preds[['c']]$actual)

ggplot(u.tdat,aes(x=c,y=pred,color=p)) + geom_point()

ggplot(u.tdat,aes(x=p,y=pred,color=c)) + geom_point()

base=ggplot(u.tdat,aes(x=c,y=p)) 
  base + geom_raster(aes(fill=pred),interpolate=TRUE)
  
  base + stat_density(aes(fill=..count..), geom='raster', position='identity')

#p.preds = data.frame(id=1:20,pred=pred[['p']],actual=preds[['p']]$actual)
#c.preds = data.frame(id=1:39,pred=pred[['c']],actual=preds[['c']]$actual)


calcplt=ggplot(preds[[d]],
               aes_string(y='est',x='id') ) + 
  geom_point() + 
  geom_errorbar(aes(ymax=up,ymin=down)) +
  geom_line(aes(y=actual))


#conditional means (observed)

u.tdat=unique(tdat[,c('a','p','c')])

tt=data.frame(p=window(u.tdat$p,1),c=window(u.tdat$c,1))
tt.prop=prop.table(table(tt),1)
View(tt.prop)
 

#@@@@@@posterior predictive and bayesian p-values

#simulate y distributions from estimates

ytilde = xhat %*% t(as.matrix(allmods[[testmodel]]$betas))

ytilde = sapply(1:1000,function(i)
              ytilde[,i] + rnorm(10000,mean=0,sd=allmods[[testmodel]]$sigma[i]))

#omnibus bayesian pvalues
print(summary(tdat$y2))
print(summary(as.vector(ytilde)))

#mean
sum(apply(ytilde,2,mean)<mean(tdat$y2))/ncol(ytilde)
#sum(apply(ytilde,2,max)<max(tdat$y2))/ncol(ytilde)
#sum(apply(ytilde,2,min)>min(tdat$y2))/ncol(ytilde)


#bayesian p-values by period 
library(dplyr)

actual.period = aggregate(tdat$y2,by=list(tdat$p),mean)
ytilde.period =  sapply(1:ncol(ytilde),function(i)
                        aggregate(ytilde[,i],by=list(tdat$p),mean)[[2]]
                        )
sapply(seq_along(actual.period[[1]]),function(i)
              sum(ytilde.period[i,]>actual.period[i,2])/nrow(ytilde)
)

plot(actual.period)

#bayesian p-values by cohort 
actual.cohort = aggregate(tdat$y2,by=list(tdat$c),mean)
ytilde.cohort =  sapply(1:ncol(ytilde),function(i)
  aggregate(ytilde[,i],by=list(tdat$c),mean)[[2]]
)
sapply(seq_along(actual.cohort[[1]]),function(i)
  sum(ytilde.cohort[i,]>actual.cohort[i,2])/nrow(ytilde)
)

plot(actual.cohort)

