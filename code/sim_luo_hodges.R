#Dev R 3.3.0 "Supposedly Educational"
#simulates data for 4 scenarios from Luo & Hodges SMR 2015; p. 14
#
#Bryce Bartlett

#clear cache
rm(list=ls())
source('config~.R')

#print(
#  round(runif(1)*1000000)
#)

set.seed=421208


#specific source functions (will become package)
#source('funs.R')

#@@helper funcitons

#generate combinations of 20 ages and 20 periods
a=rep(1:20,20)
p=a[order(a)]
c=p-a

combos=data.frame(
  a=a,a2=a^2,
  p=p,p2=p^2,c=c,c2=c^2
)

#generate "x" data
xcov=cov(combos[c('a','p','c')])
Omega=cov2cor(xcov)
combos$x.ind=
  round(runif(nrow(combos),max=max(c),min=min(c)))
  combos$x.ind2=combos$x.ind^2
#add a tiny bit of noise that breaks linear dependence
combos$x.cor = c + runif(nrow(combos),max=1,min=-1)
combos$x.cor2=combos$x.cor^2

 #rm(a,p,c)

#repeat data for 25 individuals
dat=combos; for(i in 2:25){dat=rbind(dat,combos)}; #rm(combos,i)
#confirm only 400 unique combos
cat('400 Unique Observations Test:',
  nrow(unique(dat)) == 400, '\n\n'
)

n=nrow(dat)

#id true betas for 4 scenarios (eq. 8-11)

t.beta=data.frame(
  a=0.3,a2=-0.01,p=-0.04,p2=0.02,c=0.35,c2=-0.0015
)
t.beta[2,] = t.beta[1,]; t.beta[2,c('a','a2')]=0
t.beta[3,] = t.beta[1,]; t.beta[3,c('p','p2')]=0
t.beta[4,] = t.beta[1,]; t.beta[4,c('c','c2')]=0

print(t.beta)

#draw unique y for each scenario y to y4 in a list of 100 times 
#/innefficeint memory -- should only store ys
dats=list(); inddats=list(); cordats=list()

for(N in 1:100){
  dats[[N]] = dat[,c('a','a2','p','p2','c','c2')]
  inddats[[N]] = dat[,c('a','a2','p','p2','x.ind','x.ind2')]
  cordats[[N]] = dat[,c('a','a2','p','p2','x.cor','x.cor2')]
  for(s in 1:nrow(t.beta)){
    b = as.matrix(t.beta[s,])
    x=as.matrix(dat[,1:length(b)])
    x.ind=as.matrix(inddats[[N]][,1:length(b)])
    x.cor=as.matrix(cordats[[N]][,1:length(b)])
    
   dats[[N]][,paste0('y',s)]=rnorm(n,mean=x %*% t(b),sd=1)
   inddats[[N]][,paste0('y',s)]=rnorm(n,mean=x.ind %*% t(b),sd=1)
   cordats[[N]][,paste0('y',s)]=rnorm(n,mean=x.cor %*% t(b),sd=1)
  }
}

######
#Export data for fitted plots of scenarios based on 
#true betas
######

#pull unique data -- should make this a loop
uns = list(
  c = unique(dats[[1]]$c), a=unique(dats[[1]]$a), p=unique(dats[[1]]$p)
)

pltdat=list()

#need to hold means for margins (i.e. margins command in stata)??
dims=c('a','p','c')

for(d in dims){
  o.dims = dims[which(!dims==d)]
  o.means=apply(dat[,o.dims],2,mean)
  xdat=data.frame(uns[[d]],o.means[1],o.means[2]); colnames(xdat)=c(d,o.dims)
  xdat$a2 = xdat$a^2; xdat$p2 = xdat$p^2; xdat$c2 = xdat$c^2
  x = as.matrix(xdat)
  b = as.matrix(t.beta)[,colnames(x)]
  
  #test for matching
  if(!all(colnames(x) == colnames(b))){stop('Problem with matrix ordering.')}

  pltdat[[d]] = as.data.frame(x %*% t(b)); 
  colnames(pltdat[[d]]) = paste0('s',1:4)
  
  pltdat[[d]]$id = uns[[d]] 

}

save(pltdat,file=paste0(datdir,'luo_sim_fits.RData'))


  
#####
#Save and export simulated data
#####

save(dats,file=paste0(datdir,'luo_sim.RData'))
save(inddats,file=paste0(datdir,'luo_sim_indx.RData'))
save(cordats,file=paste0(datdir,'luo_sim_corrx.RData'))

save(t.beta,file=paste0(datdir,'truebeta.RData'))

tdat = dats[[sample(seq_along(dats),1)]][,c('a','p','y1','y2','y3','y4')]

#save 1 sset of testdata
save(tdat,file=paste0(datdir,'testdat.RData'))

#write .csv of 1 set for pasting to stata
write.csv(tdat[,c('a','p','y1','y2','y3','y4')],
          "H:/projects/apc/output/sim_data/tdat.csv")

