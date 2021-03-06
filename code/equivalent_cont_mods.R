###equivalents
#written on paper reducing secenario 1 to a and p only

rm(list=ls())

#source('sim_r2.R',echo=TRUE)

age.x=1:20
p=10.5 #(held at mean for interaction

plot(age.x,-0.05*age.x-0.0115*age.x^2-0.003*age.x*p,
     type='l',lty=2,ylim=c(-7,4))
lines(age.x,.3*age.x-0.01*age.x^2)

##
#load(paste0(datdir,'sim2_tbeta.RData'))

t.beta = data.frame(a=-0.81,a2=-0.006,
                    p=0,p2=0,
                    c=-0.99,c2=-0.09)


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

#r2=0.3
#evar=(1-r2)/r2
#e=rnorm(n,0,sqrt(evar))
e=rnorm(n,0,1)

#simulate and save data
dat$s1 = as.vector(as.matrix(dat)%*%t(as.matrix(t.beta)))
dat$y1 = dat$s1 + e

#need to expand and proove
tform = function(betas,dim.from){
  #helper funciton transforms apc values from vector of t.beta (named apc)
  #input is named vector betas (a = age, p = period, c = cohort; a2=a^2, etc).
  # dim.from is a character vector 'a' or 'p' or 'c' for pulling "out" the effect
  #output is vecor with coefficients recalculated
  
  #a = p - c
  #p = a + c
  #c = p - a
  
  betas = as.matrix(betas)
  
  #select effects to transform
  fr = betas[,grepl(dim.from,colnames(betas))]
  #colnames(fr) = grepl(dim.from,colnames(betas))
  quad = grepl('2',names(fr))
  lin = !quad
  
  #linear transform
  tr =  betas[1,c('a','c')] + fr[lin]
  
  #quadratic transform
  tr2 = c(betas[1,c('a2','c2')],fr[quad]) + fr[quad]
  names(tr2) = c('a2','c2','ac') 
  
  tf = c(tr,tr2)
  leftovers = c('a','a2','c','c2')
  tf[leftovers] = tf[leftovers] + betas[,leftovers]
  
  return(as.data.frame(t(tf)/2))
  
  #polynomial trans
  #https://math.stackexchange.com/questions/454007/matrix-representation-of-polynomial
  #higher order
  #https://en.wikipedia.org/wiki/Tensor
  
}

#test transform function -- by hand
testb = data.frame(a=1,p=1,c=1,a2=1,p2=1,c2=1)
tf.testb = tform(testb,dim.from='a')

print(t.beta)
tf.beta = tform(betas=t.beta,dim.from='a')
print(tf.beta)



#direct test
dat$r.yhat = as.matrix(dat[,colnames(testb)]) %*% t(as.matrix(testb))
dat$t.yhat = as.matrix(dat[,colnames(tf.testb)]) %*% t(as.matrix(tf.testb))

#lm test
dat$y1 = dat$r.yhat + rnorm(nrow(dat),mean=0,sd=2)
summary(lm(y1~a+c+a2+c2+ac,data=dat))

#test sim_2R data
#tdat$c = tdat$p-tdat$a
#tdat = tdat %>% 
#  mutate(a2=a^2,c2=c^2,ac=a*c)

#print(summary(lm(y1~a+c+a2+c2+ac-1,data=tdat)))
#print(tform(betas=t.beta,dim.from='p'))
