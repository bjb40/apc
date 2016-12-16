#R 2.14.0
# Also for Ind Study

#@@@@@@@@@
#Globals
#@@@@@@@@@

start = Sys.time()

outdir = "../output"
subref = paste(outdir,"/subdat.csv",sep='')
cleanref = paste(outdir,"/cleandat.csv",sep='')

#@@@@@@@@@@@@ 
#Load data
#Limit to happy 2006 GSS
#@@@@@@@@@@@@

if (exists("clean") == FALSE & file.exists(cleanref)){
  clean = read.csv(cleanref,na.strings='.')
}

sub = subset(clean, year == 2006, select=c(happy,sex,recession.cohort,satfin,babies,preteen,teens))
sub [,2] = sub[,2] - 1 #recode for females 0,1

xnames = c('female','rec_cohort','satfin','babies','preteen','teens')
yname = 'happy'

#listwise delete missing 
dat = na.omit(matrix(as.numeric(unlist(sub)),nrow=nrow(sub)))
ndel = nrow(sub) - nrow(dat)
rm(sub, clean)

intercept =matrix(1,nrow(dat))
y = (dat[,1])
x = cbind(intercept,(dat[,2:ncol(dat)]))

#@@@@@@@@@@@@
#Gibbs algorithm (p. 171 Lynch)
#@@@@@@@@@@@@

# establish vectors and quantities
iter = 2000

s2=matrix(1,iter)
b= matrix(0,iter,ncol(x))

xtxi = solve(t(x)%*%x)
pars=coefficients(lm(y~x-1))

#simulate beta from mvn
for (i in 2:iter){
    
  #simulate beta from multivariate normal conditional
  b[i,]=pars+t(rnorm(length(pars),mean=0,sd=1))%*%chol(s2[i-1]*xtxi)
    
  #simulate sigma from its inverse gamma distribution
  s2[i] = 1/rgamma(1,nrow(x)/2,.5*t(y-x%*%(b[i,]))%*%(y-x%*%(b[i,])))
}

#generate posterior predictive
burn = 500
ppd = matrix(0,length(y),iter-burn)


for( i in (burn+1):iter){
  for (j in 1:nrow(ppd)){
    ppd[j,i-(burn+1)] = t(x[j,]) %*% b[i,] + rnorm(1,mean=0,sd=sqrt(s2[i]))
  }
}

post = cbind(b[(burn+1):iter,],s2[(burn+1):iter])
colnames(post) = c('int',xnames,'s2')
write.csv(post,paste(outdir,'/linear_post.csv',sep=''))

dt = cbind(y,x[,2:ncol(x)])
colnames(dt) = c(yname,xnames)
colnames(ppd) = c(1:ncol(ppd))

write.csv(cbind(dt,ppd),paste(outdir,'/linear_ppd.csv',sep=''))

post = read.csv(paste(outdir,'/linear_post.csv',sep=''))
ppd = read.csv(paste(outdir,'/linear_ppd.csv',sep=''))



#summarize posterior
apply(post[-1],2,mean)
apply(post[-1],2,quantile, probs=c(0.025,0.975))

#rsquare

