#runs 5 simulations

#####check for and create simulations folder
source('config~.R')
if(!dir.exists(paste0(outdir,'simdata/'))){
  dir.create(paste0(outdir,'simdata/'))
}

######you can add a check here to 'add' to simulations
#instead of replace ... also generate an outfile...

for(drawsimulations in 1:5){
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
  }
  
  draw_betas()   
  source('sim_r2.R') 
  #source('algorithm.R')
  library(apcwin)
  #apcsamp()
  tdat$c = tdat$p-tdat$a
  smp = apcsamp(tdat,dv='y1',method='ml',
                samples=2500,cores=4)
  eff = draw_effs(smp,tol=.001) #round @2 -- parallelize to go further
  save(list = c('eff','t.beta','pltdat'),
       file=paste0(outdir,'simdata/sim',drawsimulations,'a.RData'))
  
  
  
  #kill random dimension of effects and rerun --- all but last; leave last for example
  k = sample(c('a','p','c'),1)
  t.beta[,c(k,paste0(k,'2'))]= 0
  save(t.beta,file=paste0(datdir,'sim2_tbeta.RData'))
  source('sim_r2.R')
  tdat$c = tdat$p-tdat$a
  smp = apcsamp(tdat,dv='y1',method='ml',
                samples=2500,cores=4)
  eff = draw_effs(smp,tol=.001) #round @2 -- parallelize to go further
  save(list = c('eff','t.beta','pltdat'),
       file=paste0(outdir,'simdata/sim',drawsimulations,'b.RData'))
  
}

