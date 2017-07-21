source('config~.R')
library(ggplot2)
library(reshape2)

draw = data.frame(rdirichlet(1000,c(2,1,.5)))

colnames(draw)=c('a','p','c')
plt=ggplot(draw,aes(x=a,y=p,color=c)) + geom_point()
print(plt)

m=39

len=m

alphs.e = rep(1*m,len)
#alphs = runif(len)*alphs.e
alphs = c(1,1,1,runif((len-3),0,(len+10)))
n.w = len

dr2 = rdirichlet(10000,alphs); 
dr2 = data.frame(t(apply(dr2,1,function(x) round(cumsum(x)*n.w))))
colnames(dr2) = seq_along(alphs)
#dr2$grp = 'Alpha is 0.2, 1, or 2'
dr3 = rdirichlet(10000,alphs.e); 
dr3 = data.frame(t(apply(dr3,1,function(x) round(cumsum(x)*n.w))))
colnames(dr3) = seq_along(alphs)
#dr3$grp = 'Alpha = 1'

#```should append a 0 to the lwft instead of assign the 1...

dr2.b = dr3.b = dr3
dr2.b[,len] = dr2.b[,len] = len
dr3.b[,len] = dr3.b[,len] = len

for(c in 1:len){
  
  if(c==1){
    dr3.b[,c]=TRUE
    dr2.b[,c]=TRUE
      } else{
    dr3.b[,c]=dr3[,c] != dr3[,c-1] 
    dr2.b[,c]=dr2[,c] != dr2[,c-1] 
    }
}

dr3.row= rowSums(dr3.b[,1:len])
dr2.row = rowSums(dr2.b[,1:len])

cat('Window Breaks equal,unequal:', mean(dr3.row), mean(dr2.row))

#print('Mean probabilities.')
dr3=round(colMeans(dr2.b[,1:len]),2)
dr2=round(colMeans(dr3.b[,1:len]),2)

dt = melt(rbind(dr2,dr3)) %>%
  mutate(Var1=ifelse(Var1=='dr2','\u03b1 = [1,1,1,1,1,1,1,1,1,1]',
                '\u03b1 = [0.2,1,2,0.2,1,2,1,1,3,1]'))
dt$Var2 = factor(dt$Var2)

panel.plt = ggplot(dt, aes(y=value,x=Var2)) + 
  geom_bar(stat='identity') + 
  facet_grid(Var1~.) +
  theme_classic() + 
  xlab('') + ylab('Distribution of Probablities for a Window Break') +
  ylim(0,1)

print(panel.plt)
ggsave(paste0(imdir,'dirichelet_eg.pdf'))
