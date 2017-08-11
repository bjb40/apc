#analyzes allmods object produced by prep_gss.R
#Bryce Bartlett

rm(list=ls())

source('config~.R')

#load fit object -- should be chains
load(paste0(outdir,'empirical_res.RData'))

#reconform dat (this happens "inside" chains, so is lost) 
dat = dat %>%
  rename(a=age,
         p=year,
         c=birthyear)

#combine chains --- can calculate R-hats here...
extract = function(l,name,as.df=FALSE,span=NULL){
  #extracts and combines name object from list 'l'
  #span limits the extraction, if, for example, I need to drop first off
  if(is.null(span) & !as.df){span=1:length(l[[1]][[name]])}
  
  res=do.call(base::c,lapply(l,function(x) x[[name]][span]))
  if(as.df){
    res=do.call(rbind,lapply(l,function(x) x[[name]]))
  }
  
  return(res)
}


#post-processing of models previously run

#@@
#extracta data from all chains 
win = extract(chains,'win',as.df=TRUE)  
allmods = extract(chains,'allmods',span=2:250)
xhats = extract(chains,'xhats',span=2:250)
effects = extract(chains,'effects',span=2:250)
n.samples = sum(unlist(extract(chains,'n.samples',as.df=TRUE)))
acc = sum(unlist(extract(chains,'acc',as.df=TRUE)))
bound = sum(unlist(extract(chains,'bound',as.df=TRUE)))

#breaks have a wierd pattern -- extracted differently
breaks = list(a=NULL,p=NULL,c=NULL)
for(l in names(breaks)){
  breaks[[l]] = do.call(base::c,lapply(chains,function(x) x[['breaks']][[l]][2:250]))
}

#@@
#remove chain object
rm(chains)

bics=unlist(lapply(allmods,FUN=function(x) x$bic))
bics_prime=unlist(lapply(allmods,FUN=function(x) x$bic_prime))
r2=unlist(lapply(allmods,FUN=function(x) mean(x$r2)))
  
best=which(bics_prime==min(bics_prime))

  
  library(ggplot2)
  library(gridExtra)
  
  best.plt = list()
  preds = list()
  
  for(d in c('a','p','c')){
    
    preds[[d]]=as.data.frame(apply(effects[[best]][[d]],2,mean))
    colnames(preds[[d]])='est'
    rng=apply(effects[[best]][[d]],2,quantile,c(0.025,0.975))
    preds[[d]]$up = rng[2,]
    preds[[d]]$down = rng[1,]
    #s1-s4 are for scenarios --- needs to match with y1-y4
    
    #pldat doesn't exist -- need to figure about ordering
    preds[[d]]$id=unique(dat[,d])[order(unique(dat[,d]))]

    
  }
  
  
  ##post-processing -- model averaging
  #averaging algorithm from ... eq 35 from Rafferty SMR
  #http://scicomp.stackexchange.com/questions/1122/how-to-add-large-exponential-terms-reliably-without-overflow-errors
  #http://stats.stackexchange.com/questions/249888/use-bic-or-aic-as-approximation-for-bayesian-model-averaging
  #fixing overflow issue from suggestion above after logging to take difference
  
  #w = exp(-.5*bics)/sum(exp(-.5*bics)))
  #w_prime=exp(-.5*bics_prime)/sum(exp(-.5*bics_prime))
  
  #sbics = bics[order(bics)][1:10]
  
  k=min(bics)
  d=-.5*(bics-k)
  w=exp(d)/sum(exp(d))
  
  k = min(bics_prime)
  d=-.5*(bics_prime-k)
  w_prime=exp(d)/sum(exp(d))
  
  
  #add weight to apc windows dataframe
  win$wt=w; win$w_prime=w_prime; win$r2=r2; win$bic=bics; win$bic_prime=bics_prime
  win$rmse = unlist(lapply(allmods,FUN=function(x) mean(x$rmse))) 
  #inverse weight by rmse (larger values are smaller weights)
  win$rmsewt = 1/win$rmse/sum(1/win$rmse)
  win$mse = unlist(lapply(allmods,FUN=function(x) mean(x$rmse^2))) 
  win$msewt=1/win$mse/sum(1/win$mse)
  
  win$modnum=1:nrow(win)
  
  effects = effects[2:length(effects)]
  
  #select weight
  #use.wt='rmsewt'
  use.wt='wt'
  
  for(mod in seq_along(effects)){
    effects[[mod]]$w=win[,use.wt][mod]
  }
  
  #print weighted mean of windows...
  print(
    apply(win[,c('a','p','c')],2,weighted.mean,w=win[,use.wt])
  )
  
  #weighted mean
  for(l in seq_along(effects)){
    for(t in c('a','p','c')){
      if(is.null(effects[[l]][[t]])){
        effects[[l]][[t]] = 0
      }
    }
  }
  
  wtmn=function(l,var,w){
    #makes weighted mean
    #l is list results
    #var is charactername
    r=NA
    if(length(l[[var]])==1){r=0} else{
      r=apply(l[[var]],2,mean)*w
    }
    return(r)
  }
  
  wtquant=function(l,var,w,q){
    #takes weighted mean of quantile
    #l is list results
    #var is charactername
    #q is quantile
    r=NA
    if(length(l[[var]])==1){r=0} else{
      r=apply(l[[var]],2,quantile,q)*w
    }
    return(r)
  }
  
  
  mean.plt=list()
  
  for(d in names(preds)){
    
    preds[[d]]$m_est =  rowSums(
      as.data.frame(
        lapply(effects,FUN=function(e)
          wtmn(e,d,e$w))
      ))
    
    
    
    preds[[d]]$m_down=rowSums(
      as.data.frame(lapply(
        effects,FUN=function(e)
          wtquant(e,d,e$w,0.025))
      )
    )
    
    preds[[d]]$m_up=rowSums(
      as.data.frame(lapply(
        effects,FUN=function(e)
          wtquant(e,d,e$w,0.975))
      )
    )
    
    
    
  }
  
  
  df = do.call(rbind,preds)
  df$type=substr(rownames(df),1,1)
  
  b.plt=ggplot(df,aes(x=id,y=est)) +
    geom_line(alpha=0.5) +
    geom_ribbon(aes(ymin=down,ymax=up), alpha=0.5)  +
    facet_grid(.~type,scales='free_x') +
    theme_minimal() +
    ylab('Effect') +
    xlab('APC Value')
  
  print(b.plt)
  
  m.plt=ggplot(df,aes(x=id,y=m_est)) +
    geom_line() +
    geom_ribbon(aes(ymin=m_down,ymax=m_up), alpha=0.25) +
    facet_grid(.~type,scales='free_x') +
    theme_minimal() +
    ylab('Effect') +
    xlab('APC Value')
  
  print(m.plt)
  
  #save entire environment for post-processing
  