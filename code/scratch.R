########
#for new dgp


est = apcsamp(dat=tdat,dv='y',
              cores=4,method='ml',
              chains=4,samples=25)

summary(est)

effs = draw_effs(est)



###
#real data

#pool real effects
r.effs = list(a=a.eff,
              p=pdraw$var.eff,
              c=cdraw$var.eff)

#pull unique values for plotting
r.effs=lapply(names(r.effs),function(x){
  df = data.frame(eff=r.effs[[x]],val=tdat[,x])
  df$dim=x
  return(unique(df))
})  

r.effs = do.call(rbind,r.effs)

plot(effs) +
  geom_line(data=r.effs,
            aes(x=val,y=eff),linetype=2) + 
  theme_classic()

ggplot(r.effs,aes(x=val,y=eff)) + 
  geom_line(linetype=2) + 
  theme_classic() +
  facet_wrap(~type,scales='free_x')

