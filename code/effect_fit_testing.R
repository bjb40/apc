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

#conditional means (observed)

u.tdat=unique(tdat[,c('a','p','c')])

tt=data.frame(p=window(u.tdat$p,1),c=window(u.tdat$c,1))
tt.prop=prop.table(table(tt),1)
View(tt.prop)



