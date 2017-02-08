#omitted variable tests

library(dplyr)

tmod = list()

tmod[[1]] = lm(y2~p+c+I(p^2) + I(c^2),tdat)

p.mar = tdat %>% group_by(p) %>% summarize(c=mean(c)) %>% select(p,c) %>% ungroup
p.pred = predict.lm(tmod[[1]],p.mar)

ggplot() + geom_line()


predict.lm(tmod[[1]],data.frame(p=1:20,c=rep(0,20)))

ndat= unique(tdat[,c('p','c')])
View(predict.lm(tmod[[1]],ndat))


#conditional bivariate
#https://onlinecourses.science.psu.edu/stat414/node/118
#note that you may just calculate sets of different functional forms