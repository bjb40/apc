
tdat = simulated$simdat

tt = do.call(rbind,
       lapply(sampobj$breaks[['p']], function(x)
       as.numeric(unique(tdat$p) %in% x))
)
colnames(tt) = unique(tdat$p)
lp = length(unique(tdat$p))

#tt[1,]
#sampobj$breaks[['p']][[1]]
chainlen = nrow(tt)/sampobj$chains
chain = factor(ceiling(1:nrow(tt)/chainlen))
tt =split(as.data.frame(tt),chain)

tt = lapply(tt,mcmc)
tt = as.mcmc.list(tt)

gelman.plot(tt[,6, drop=FALSE])
gelman.plot(tt[,7, drop=FALSE])


gelman.diag(tt[,1:(lp-1),drop=FALSE])
autocorr.diag(tt)
