#runs 5 simulations
#identifies single window breaks
#bryce bartlett

rm(list=ls())

#####check for and create simulations folder
source('config~.R')
if(!dir.exists(paste0(outdir,'stepsim/'))){
  dir.create(paste0(outdir,'stepsim/'))
}

######you can add a check here to 'add' to simulations
#instead of replace ... also generate an outfile...

for(drawsimulations in 1:1){
  cat('Simulation',drawsimulations,'...\n')
  #funciton is nested because sources below delete all in memory
  source('config~.R')
  draw_betas = function(){
    #draws random set of betas
    t.beta=data.frame(t(runif(6,-1,1)))
    colnames(t.beta)=c('a','a2','p','p2','c','c2')
    #rescale quadratic effects to be smaller
    t.beta[,c('a2','p2','c2')] = t.beta[,c('a2','p2','c2')]/10
    #print(t.beta)
    save(t.beta,file=paste0(datdir,'sim2_tbeta.RData'))
    return(t.beta)
  }
  
  draw_betas()   
  source('sim_r2.R') 
  #source('algorithm.R')
  library(apcwin)
  #apcsamp()
  tdat$c = tdat$p-tdat$a
  
  #####
  #prepare a loop that will march through all combinations
  dims = c('a','p','c')
  
  #set a list that can be used to make window breaks
  rng = lapply(dims,
               function(x){ 
                 r=range(tdat[,x])
                 return((r[1]-1):r[2])
                 })

  names(rng) = dims
    
  #generate unique combinations of window breaks marching
  #one by one by one ...
  
  delone = function(lst){
    #helper function to return a list of window breaks with
    #one deleted
    res=list()
    #note that the first and last window breaks are undefined 
    for(i in 2:(length(lst)-1)){
      del = lst[i]
      res[[i-1]] = unlist(lst[!lst==del])
    }
    
    return(res)
    
  }#end helper function
  
#generate list of unique window breaks for each dimension
  
brks = lapply(rng,delone)

#cycle through and create a dataframe for each unique combo

#make frame for "true data"
x = tdat

###
#create objects to use the draw_effects method in apcwin
n.samples = 0
win = data.frame(a=as.numeric(NA),
                 p=as.numeric(NA),
                 c=as.numeric(NA))

modsum = data.frame( r2=as.numeric(NA),
                     sigma=as.numeric(NA),
                     bic=as.numeric(NA),
                     bic_prime=as.numeric(NA))

breaks=list(a=NULL,p=NULL,c=NULL)

#total sample space
tots = prod(unlist(lapply(brks,length)))

for(i.a in seq_along(brks[['a']])){
  for(i.p in seq_along(brks[['p']])){
    for(i.c in seq_along(brks[['c']])){
      n.samples=n.samples+1
      
      breaks[['a']][[n.samples]] = brks[['a']][[i.a]]
      breaks[['p']][[n.samples]] = brks[['p']][[i.a]]
      breaks[['c']][[n.samples]] = brks[['c']][[i.a]]
      
      
      x$a = window(tdat$a,breaks=brks[['a']][[i.a]])
      x$p = window(tdat$p,breaks=brks[['p']][[i.p]])
      x$c = window(tdat$c,breaks=brks[['c']][[i.c]])
  
  #run models and store
  xmat = model.matrix(~a+p+c,data=x)
  y=x$y1
  m = lin_ml(y=y,x=xmat)
 
  modsum = rbind(modsum,
                 m[c('sigma','r2','bic','bic_prime')])
  nr=data.frame(a=length(levels(x$a)),
                p=length(levels(x$p)),
                c=length(levels(x$c)))
  
  
  win=rbind(win,nr)

#from chain sampler --- (so you can use effects method)   
#  res = list(
#    modsum=modsum[2:nrow(modsum),],
#    breaks=breaks,
#    win=win[2:nrow(modsum),],
#    n.samples=n.samples,
#    acc=acc,
#    bound=bound,
#    method=method
#  )
  
  #mods[[nm]] = m
#  if(n.samples%%100 ==0){
#    print(paste('Sampled',n.samples,'of',tots))}
  
    }#end unique c combos
  }#end unique p combos
}#end unique a combos

  

modsum=modsum[2:nrow(modsum),]

k=min(modsum$bic)
d=-.5*(modsum$bic-k)
modsum$w=exp(d)/sum(exp(d))

k = min(modsum$bic_prime)
d=-.5*(modsum$bic_prime-k)
modsum$w_prime=exp(d)/sum(exp(d))

stepsamp = list(
      data = tdat,
      modsum=modsum,
      breaks=breaks,
      summaries=win[2:nrow(win),],
      n.samples=n.samples-1,
      acc=1,
      bound=0,
      method='ml',
      dv='y1',
      apc=c('a','p','c')
    )

class(stepsamp) = append(class(stepsamp),'apcsamp')

#summary(stepsamp)

eff = draw_effs(stepsamp,tol=0.001)



#save.image(paste0(outdir,'steve_algo_sim3.RData'))

save(list = c('eff','t.beta','pltdat'),
     file=paste0(outdir,'stepsim/sim',drawsimulations,'.RData'))

} # end simulator

