#runs 15 simulations

draw_betas = function(){
#draws random set of betas
  t.beta=data.frame(t(runif(6,-1,1)))
  colnames(t.beta)=c('a','a2','p','p2','c','c2')
  #rescale quadratic effects to be smaller
  t.beta[,c('a2','p2','c2')] = t.beta[,c('a2','p2','c2')]/10
  #print(t.beta)
  save(t.beta,file=paste0(datdir,'sim2_tbeta.RData'))
}

#for(drawsimulations in 1:15){
for(drawsimulations in 1:3){
  draw_betas()
  source('sim_r2.R')
  source('algorithm.R')
  
  #kill random dimension of effects and rerun
  k = sample(c('a','p','c'),1)
  t.beta[,c(k,paste0(k,'2'))]= 0

}

