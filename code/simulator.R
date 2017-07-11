#runs 5 simulations


for(drawsimulations in 1:5){
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
  source('algorithm.R')
  
  #kill random dimension of effects and rerun --- all but last; leave last for example
  k = sample(c('a','p','c'),1)
  t.beta[,c(k,paste0(k,'2'))]= 0
  save(t.beta,file=paste0(datdir,'sim2_tbeta.RData'))
  source('sim_r2.R')
  source('algorithm.R')
}

