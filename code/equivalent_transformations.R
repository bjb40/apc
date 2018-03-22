
t.beta = data.frame(a=1,a2=1,
                    p=0,p2=0,
                    c=1,c2=1)

##########
#simulate data
a=rep(1:20,20)
p=a[order(a)]
c=p-a

combos=data.frame(
  a=a,a2=a^2,
  p=p,p2=p^2,c=c,c2=c^2
)

dat=combos; for(i in 2:25){dat=rbind(dat,combos)}; #rm(combos,i)
#confirm only 400 unique combos
cat('400 Unique Observations Test:',
    nrow(unique(dat)) == 400, '\n\n'
)

n=nrow(dat)

#r2=0.003
#evar=(1-r2)/r2
#e=rnorm(n,0,sqrt(evar))
e=rnorm(n,0,1) # rmse of 100 gives r-sq of 0.7; doesn't work as well

#simulate and save data
dat$s1 = as.vector(as.matrix(dat)%*%t(as.matrix(t.beta)))
dat$y1 = dat$s1 + e

#####
#test
summary(lm(y1~a+I(a^2)+c+I(c^2),data=dat))

########
#algebraic transformation
#  a + a^2 + c + c^2 
#= (p-c) + (p-c)^2 + c^2
#= p - c + p^2 - 2pc + c + c^2 + c^2
#= p + p^2 -2pc + 2c^2

summary(lm(y1~p + I(p*c) + I(p^2) + I(c^2),data=dat))

