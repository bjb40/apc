##42 is single windows for a and p with no c y4 

#apc_lm(y4~a+p,tdat,age='a',per='p')

for(d in c('a','p')){
  
  preds[[d]]=as.data.frame(apply(effects[[42]][[d]],2,mean))
  colnames(preds[[d]])='est'
  rng=apply(effects[[42]][[d]],2,quantile,c(0.025,0.975))
  preds[[d]]$up = rng[2,]
  preds[[d]]$down = rng[1,]
  #s1-s4 are for scenarios --- needs to match with y1-y4
  #preds[[d]]$actual=pltdat[[d]]$s1[order(pltdat[[d]]$id)]
  preds[[d]]$actual=pltdat[[d]]$s4[order(pltdat[[d]]$id)]
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


grid.arrange(best.plt[['a']],
             best.plt[['p']],
#             best.plt[['c']],
             ncol=3)

#x estimates
x = tdat[,c('a','p')]
x$a=window(x$a,winlength=1)
x$p=window(x$p,winlength=1)

xmat=model.matrix(~.,data=as.data.frame(x))
xhat = xhats[[42]]

summary(lm(tdat$y4~xmat))

rev=list(a=apply(round(xhat[['a']]) %*% t(as.matrix(allmods[[42]]$betas)),1,mean),
         p=apply(round(xhat[['p']]) %*% t(as.matrix(allmods[[42]]$betas)),1,mean))

plot(1:20,rev[['a']]); lines(1:20,preds[['a']]$actual)

plot(1:20,rev[['p']]); lines(1:20,preds[['p']]$actual)


