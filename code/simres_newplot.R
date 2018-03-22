
require(apcwin)

dts = list.files('H:/projects/apc/output/simdata/',
                 full.names=TRUE)

ftplots = list()
apcdim = list()
id_dim = list()

for(f in seq_along(dts)){
  
  load(dts[[f]])

  
######
#you can loop to load here...

real.dat = lapply(names(pltdat),
               function(nm){
                 x = as.data.frame(pltdat[[nm]])
                 x$x = as.numeric(row.names(x))
                 x$dim = nm
                 return(x)
               })

real.dat = do.call(rbind,real.dat)

ftplots[[f]] = plot(eff) + 
  geom_line(data=real.dat,aes(x=id,y=s1),
            lty=2) +
  facet_wrap(~dim,scales='free_x') +
  theme_void()

#calculate overlap with 95% CI
ld = ftplots[[f]]$data %>%
  rename(id=x)

#NAMES FOR LL & UL ARE BACKWARDS -- NEED TO FIX IN PACKAGE!!
ld = merge(ld,real.dat,by=c('id','dim'))
ld = ld %>%
  mutate(ol = ifelse(s1<ll & s1>ul,1,0))

summary(eff$sampobj)
cat('\n\n\n\n')

post.comp = length(unique(eff$sampled))

cat('Posterior is composed of',post.comp,'models.\n\n')
#test against limit case

ll = eff$sampobj$limits
apc_bic = mean(eff$fit$bic)
rdim_bic = lapply(ll,function(x) x$bic)

abic = c(apc_bic,rdim_bic)
names(abic)[1] = 'apc'
abic = unlist(abic)
mnb = min(abic)

awt = exp(-0.5*(abic-mnb))/sum(exp(-0.5*(abic-mnb)))
print(round(awt,3))
print(t.beta)

meancol=colMeans(eff$sampobj$summaries[eff$sampled,1:3])
names(meancol) = paste0('wb.',names(meancol))
names(post.comp)='post.comp'
id_dim[[f]] = unlist(c(t.beta,awt,post.comp,mean(ld$ol),sum(meancol)))

apcdim[[f]] = all(apc_bic<rdim_bic)


if(!all(apc_bic<rdim_bic)){

  ftplots[[f]] = ftplots[[f]] +
    theme(panel.background = 
            element_rect(fill='grey',color='grey'))
    
}



#Sys.sleep(5)


}


library(gridExtra)

pdf(paste0(imdir,'new_allfits.pdf'))
  grid.arrange(grobs=ftplots)
dev.off()

View(round(rtabs,2))
grid.arrange(grobs=ftplots)


rtabs = as.data.frame(do.call(rbind,id_dim))
library(knitr)

###relation between geom_plot and posterior compostion
#suggesting multiple targets to me...
#v12 is overlap
colnames(rtabs)[12:13] = c('olap','tot.breaks')

ggplot(rtabs,aes(x=post.comp,y=olap,color=apc)) + geom_point()

#v13 is # total av window breaks
ggplot(rtabs,aes(x=tot.breaks,y=olap,color=apc)) + 
  geom_point()


sink(paste0(outdir,'sim_table.md'))
  kable(rtabs,digits=2)
sink()

#gb = lin_gibbs(ll$no_a$y,ll$no_a$x)
#ft = eff$fit$bic
#mnb = min(ft)
#ft.w = exp(-0.5*(ft-mnb))/sum(exp(-0.5*(ft-mnb)))
#print(range(ft.w)*length(ft))


################
#calculating ppd

#y = eff$sampobj$data[,eff$sampobj$dv]
#y.tilde = 

#################

###################
#drawing stuff out

allb = lapply(eff$sampobj$breaks,function(x)
  x[eff$sampled])

rng=lapply(c('a','p','c'),
           function(x) range(eff$sampobj$data[,x]))
names(rng)=c('a','p','c')

brks = list()
for(d in names(allb)){
  brks[[d]] = lapply(allb[[d]], function(x)
      rng[[d]][1]:rng[[d]][2] %in% x)
  brks[[d]] = do.call(rbind,brks[[d]])
}

#because of the way it works, the n for breaks 
#means the probability that n+1 is new effect
brks.sum = lapply(brks,colMeans)

