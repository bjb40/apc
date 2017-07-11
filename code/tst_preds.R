##testing plots

sdat = tdat
sdat$a.w = window(tdat$a,breaks=breaks$a[[best]])
sdat$p.w = window(tdat$p,breaks=breaks$p[[best]])
sdat$c.w = window(tdat$c,breaks=breaks$c[[best]])

sdat$a2 = sdat$a^2
sdat$p2 = sdat$p^2
sdat$c2 = sdat$c^2


prd = function(var,o.means){
  #var is variable (text)
  #o.means is means of other variables, named list
  
  v = unique(sdat[,var])
  
  block = data.frame(v,o.means[1],o.means[2])
  colnames(block) = c(var,names(o.means))
  block$a2 = block$a^2
  block$p2 = block$p^2
  block$c2 = block$c^2
  
  block = block[,names(t.beta)]
  
  pred = as.matrix(block) %*% t(as.matrix(t.beta))
  df = data.frame(v,pred); colnames(df) = c('id','fit')
  
  return(df)
  
}

mga = window(mean(sdat$a),breaks=attr(sdat$a.w,'breaks'))
ma=mean(sdat$a[sdat$a.w==mga])
#ma = sdat %>% group_by(a.w) %>% summarize(a=mean(a)) %>% ungroup
#ma = mean(ma$a)

mgp=window(mean(sdat$p),breaks=attr(sdat$p.w,'breaks'))
mp=mean(sdat$p[sdat$p.w==mgp])
#mp = sdat %>% group_by(p.w) %>% summarize(p=mean(p)) %>% ungroup
#mp = mean(mp$p)

#mgc = window(mean(sdat$c),breaks=attr(sdat$c.w,'breaks'))
#mc=mean(sdat$c[sdat$c.w==mgc])
mc = sdat %>% group_by(c.w) %>% summarize(c=mean(c)) %>% ungroup
mc = mean(mc$c)

new=prd('a',list(c=mc,p=mp))

print(
  best.plt$a +
    geom_line(data=new,aes(x=id,y=fit),lty=2) 
    )

"
new=prd('c',list(a=ma,p=mp))
print(best.plt$c +
        geom_line(data=new,aes(x=id,y=fit),lty=2)
      )
"
"
new=prd('p',list(a=ma,c=mc))
print(best.plt$p +
        geom_line(data=new,aes(x=id,y=fit),lty=2)
      )
"


#ttm stuff -- basis funcitons??
library(reshape2)
tt = as.data.frame(t(effects[[best]]$a))
tt$a = 1:20
ttm = melt(tt,id='a')
ttm$a2=ttm$a^2

ttlm=lm(value~a+I(a^2),ttm)
print(summary(ttlm))
print(t.beta)

xhat = data.frame(a=1:20)
xhat = cbind(predict.lm(ttlm,newdata=xhat,interval='prediction'),xhat)

print(
  best.plt$a +
        geom_line(data=new,aes(x=id,y=fit),lty=2) +
        geom_line(data=xhat,aes(x=a,y=fit),lty=3) +
        geom_ribbon(data=xhat,aes(x=a,y=fit,ymin=lwr,ymax=upr),alpha=0.25)
)


#predicted variance of beta estimator from simulated data with sim corr
#won't work, singular, duh
#sg = cov(sdat[,names(t.beta)])
#mu = colMeans(sdat[,names(t.beta)])
#n=nrow(sdat)
#k=nrow(sg)
#r_chol = t(chol(sg))
#xsim = matrix(rnorm(n,mean=0,sd=1),n,k)%*% r_chol + t(matrix(mu,k,n))


#sigma^2 xtx-1
#won't work, singular,duh
#txx = t(xsim) %*% xsim
#vcov = evar * solve(txx)

