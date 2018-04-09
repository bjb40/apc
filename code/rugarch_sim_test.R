library(rugarch)

arch1.spec = ugarchspec(variance.model = list(garchOrder=c(1,0)), 
                        mean.model = list(armaOrder=c(0,0)),
                        fixed.pars=list(mu = 0, omega=0.1, alpha1=0.8))
class(arch1.spec)
arch1.spec

set.seed(123)
arch1.sim = ugarchpath(arch1.spec, n.sim=1000)

y = attr(arch1.sim,'path')[['seriesSim']]
plot(y,type='l')