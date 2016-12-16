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
sub [,2] = sub[,2] - 1 #recode for males 0,1
#listwise delete missing

dat = na.omit(matrix(as.numeric(unlist(sub)),nrow=nrow(sub)))
ndel = nrow(sub) - nrow(dat)
rm(sub, clean)

intercept =matrix(1,nrow(dat))
y = (dat[,1])
x = cbind(intercept,(dat[,2:ncol(dat)]))

#@@@@@@@@@@@@
#Gibbs algorithm (p. 173 Lynch)
#@@@@@@@@@@@@

# establish vectors and quantities
iter = 1000

s2=matrix(1,iter)
b= matrix(0,iter,ncol(x))

xtxi = solve(t(x)%*%x)
pars=coefficients(lm(y~x-1))

#simulate sigma from inverse gamma marginal
s2 = 1/rgamma(iter,nrow(x)-ncol(x)/2,.5*t(residuals(lm(y~x-1)))%*%residuals(lm(y~x-1)))

#set ppd
ppd = matrix(0,iter,length(y))

#simulate beta from mvn
for (i in 1:iter){
  b[i,]=pars+t(rnorm(length(pars),mean=0,sd=1))%*%chol(s2[i]*xtxi)

  #generate posterior predictive
  for (j in 1:nrow(ppd)){
    ppd[j,i] = t(x[j,]) %*% b[i,] + rnorm(1,mean=0,sd=s2[i])
  }


}

#ppd checks and plots



