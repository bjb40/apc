draw = data.frame(rdirichlet(1000,c(1,100,1)))

colnames(draw)=c('a','p','c')
plt=ggplot(draw,aes(x=a,y=p,color=c)) + geom_point()
print(plt)

#alphs = c(1,1,1,1,1,5,5,5,5,.1,.1,.1,.1,.1)
alphs = c(.1,.1,.1,1,1,5,5,5,5,1,1,6,6,.1)



dr2 = data.frame(rdirichlet(1000,alphs))
dr2.sum = apply(dr2,1,function(x) floor(cumsum(x)*length(alphs)))
plot(rowMeans(dr2.sum))



