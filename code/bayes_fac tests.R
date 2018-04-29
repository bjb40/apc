

#bf=n*log(sum(residuals(m1)^2)/sum(residuals(m2)^2)) + (k1-k2)*log(n)

#' Calculates a bayes factor for an lm anova.
#' 
#' @param m1 first anova model
#' @param m2 second anova model
#' @return A number indicating the Bayes Factor
#' @references https://arxiv.org/abs/1710.02351
#' 
bf.anova = function(m1,m2){
  
  n = nrow(m1$model)
  
  if (n != nrow(m2$model)){
    stop('Different dimensions of model matrix. Confirm you are using the same data.')
  }
  
  k1 = m1$rank
  k2 = m2$rank
  sse1 = sum(residuals(m1)^2)
  sse2 = sum(residuals(m2)^2)
  
  return(n*log(sse1/sse2) + (k1-k2)*log(n))
  
}

bst = which(est$summaries$bic == min(est$summaries$bic))

tdat$pf = drop(tdat$pf)

prior = lm(y~1,data=tdat)
rlm = lm(y~a + I(a^2) + pf + cf, data=tdat)
rlmf = lm(y~factor(a) + pf + cf, data=tdat)

bf.anova(rlm,rlmf)
bf.anova(prior,rlm)

bstbics = vector()
bfs = vector()

for(i in 1:50){
print(i)
if(i>1){prior = nm}
alph = list(a = est$alphas$a[bst,],
            p = est$alphas$p[bst,],
            c = est$alphas$c[bst,])

tdat$ia = window.sample(tdat$a,alph$a)
tdat$ip = window.sample(tdat$p,alph$p)
tdat$ic = window.sample(tdat$c,alph$c)

nm = lm(y~ia+ip+ic,data=tdat)

bstbics = c(bstbics,BIC(nm))
bfs = c(bfs,bf.anova(prior,nm))

}

plot(bstbics,bfs)
which(bfs == min(bfs))
which(bfs == max(bfs))
which(bstbics == min(bstbics))

