---
author:
	- name: Bryce Bartlett
	- name: Stephen Vaisey
title: Approximating Age-Period-Cohort Estimates by Averaging Multiple Models with Varying Window Restrictions 
---

#Abstract



Key words: Computational Sociology, Methods, Life Course, Bayesian, Culture

#Introduction


#Background

Both classic and modern social theory have proposed different effects across three social dimensions of time: age, period, and cohort [@ryder_cohort_1965]. These dimensions of time operate on individuals in different ways as they flow through social space. Age effects are driven either by the common experiences associated with the biological process of aging ```citation```, or persistent age-structuring institutions, like high school or retirement ```life course citation```.  Period effects are responses of everyone to contemporaneous social experiences, like recessions or wars ```citation```. Finally, cohort effects are socialization effects. They are broad-based historical events that stick to the populations experiencing them. The great depression undoubtadly produced a period effect, but cohorts who experienced the great depression took it with them; they were qualitatively different than other cohorts for the rest of their lives. Cohorts are virtually always defined in terms of birth year.

These three dimensions of time sit at the center of some recurrent debates. Cultural sociologists are divided over whether the main driver of culture is a period process (cultural fragmentation) or a cohort process (acquired dispositions) [vaisey_cultural_2016]. ```find a couple of health things, etc```

Unfortunately, it is not trivial to estimate the unique effects of age, period, and cohort to find evidence that may settle these debates. The fundamental problem is one of identification: these three dimensions are linearly dependent. Two of the dimensions define the third. If a researcher knows an individual's age, and the year of the survey, cohort is also defined. Because these variables are exactly colinear, they are not estimable using classic statistical techniques. A number of solutions have been proposed to address this conundrum, these include: ....

The past several years have seen a resurgence in debates surrounding these issues. The thrust of the methods outlined above is to find *the* best-fitting model. In reality, modeling constraints introduce assumptions that can produce drastically different, and these assumptions are often not apparent. This is similar to situations with high multicolinearity (say in the case of income and education), where inclusion of both variables in the model often produces results that are difficult to interperet ```see if you can find something better than the stat blog on this```.

Instead of focusing on finding the *best* fitting model, as complicated as that may be, we fall back to a well-known strategy, window (or block) constraints. These models make a set of simple and straightforward assumptions that, if true, moot the identification problem. In particular, they assume that *some* of the unique ages, periods, and cohorts have identical effects. This equality blocks the age, period, or cohort together, and breaks the exact APC identity. The major criticism of these models is that a constraint in one of the variables induces unkown (and difficult to test constraints) in the other estimators []. Instead of relying on an arbitrrary set of constraints, our averages over thousands of constraints, weighted by an approximation to the probability of the model ```terminology``` using Bayesian Model Averaging (BMA) algorithms.

#Description of the Method

```revise here```

1. 
2. 
3. 

##Blocking Windows: Piecewise Constant Functions

```intro on piecewise constant functions```

In terms of window constraints, the target model, which is inestimable, is built by estimating each unique value of age, period, and cohort, as a dummy variable:

$$
E(Y) = \beta_{00} + \sum_{\lambda=2}^{\lambda} \beta_{1\lambda}A_{\lambda} +  \sum_{\rho=2}^{\rho} \beta_{2\rho}P_{\rho} +  \sum_{\kappa=2}^{\kappa} \beta_{3\kappa}C_{\kappa}
$$

where $\beta$ is an estimated effect, $\lambda$, $\rho$ and $\kappa$ index unique values for age, period, and cohort, and $A$, $P$, and $C$ stand for matricies of dummy variable series for age, period, and cohort. This model is unidentified, because, as with the continuous case, the dummy variables in any two of the matricies above condition the third matrix. In other words, the indicator variable in $C$ is a function of the indicators in $A$ and $P$ (in mathematical terms,the probability that any given cohort dummy variable is one or zero is exctly dependent on the values of A and P, so that $p(C_{\kappa}=0|A_{\lambda},P_{\kappa}$) is either exaclty zero, or exactly one).

This dependency is broken, however, by combining one or more of the unique values together so that they share the same dummy variable series. By way of example, we can construct an age dummy variable series where $A$ is sliced into two groups based on some cut-point, say $g$, so that ...

```work here```


This would identify a binary dummy variable with an older and a younger group. We can generalize this expression to an arbitrary vector of cut-points, $G$ with subscript $\gamma$ so that 


$$
\mbox{for} \gamma < max(\gamma):
A_{\lambda}  =
\begin{cases}
  1, & \mbox{if} & a > G_{\gamma} & \mbox{and} & a \leq G_{\gamma+1} \\
  0, & \mbox{if} & a < G_{\gamma} & \mbox{or} & a > G_{\gamma+1}
\end{cases}
$$


If the vector $G$ and $\lambda$ have the following properties: $max(G) = max(a)$, $min(G) < min(a)$, $\lambda \in  \{1,2, ... {\gamma-1}\}$, and $\gamma>2$, then $G$ can describe any posible sets of window restrictions for age. Generalizing cross all dimensions of APC, permuting three similar vectors (say $G^{(d)}$) will describe any model for any possible window constraints detailed in equation 1. Accordingly, all permutations of G, as defined above, constitute the model space of window constraints ($\mathscr{M}$). For any given set of data, $\matscr{M}$ is finite, but it can be quite large. For example, 10 unique ages, periods, and cohorts present ```__``` possible models. ```describe how``` Of these```__``` models, only 1 is not estimable because of perfect colinearity. This target model is (theoretially) the least biased, although it is not the most parsimonious model. The question is how to best use infromation from some subset fo possible models in $\mathscr{M}$ to estimate unbiased APC effects. Bayesian Model Averaging (BMA) provides a straightforward way to combine models. The theoretical backdrop applies to this curcumstance quite well. In principal, there is no one true (or best) model; instead estimates are conditional on models from the modeling space ($\mathscr{M}$), and have  a posterior distribution, which is calculated as a weighted average of all models [@rafferty_1995, p. 144-145].  It has been applied in diverse areas from weather forecasting to biology to social science [@fragaso_bayesian_2015] .

##Bayesian Model Averaging (BMA) and the MC3 Algorithm 

Markov Chain methods for Bayesian Model Averaging (BMA) provides a sensible way to sample over a subset of continually better-fitting models and combine their estimates to produce an approximation of APC effects. The next two sections describe the MC3 agorithm developed for BMA, and outline the unique implementation of MC3 for this particular set of models, drawing from the Dirichelet distribution. 

##Using the Dirichlet Distribution to Sample Window Groups (G)

We simulate $G$ using two sets of nuisance parameters. For each dimension ($d$) of APC, the window breaks, $G^{(d)}$ of equaiton 2 are decomposed into (1) a vector of weights on the simplex ($B$), and (2) a scalar integer, $w^d$. $G^{(d)}$ is simply the product of $B$ times $w$. We sample $w$ from the uniform distibution, as follows:

$$
w^{(d)} \sim U(2,max(\omega_{d}))
$$

Where $\omega_{d}$ is the index number for the APC effects (the $\lambda$, $\rho$, and $\kappa$ of equaiton 1). We sample the weights using the dirichelet distribution:

$$
B^{(d)} \sim Dir(\alpha_{\omega_d})
$$

Where $G$ is the matrix described above in equation 2, and the $\alpha$ for the Dirichelet distribution are for each of the vunique values in the dimension ($d$) for APC. The Dirchelet distribution draws a random projection on the Simplex, *i.e.* avector of weights $\omega_d$ weights between 0 and 1, that sum to 1 ($\sum B = 1$). Multiplying this vector by the scalar $w$ provides a set of numbers which sum to $w^{(d)}$, and represents the number of windows in that dimension. We reduce the siplex by selecting only the unique values of the cumulative sums of $B$ times the mean length of the window $w$ as follows:

$$
G^{(d)}_ = \Bigg \lfloor w^{(d)} \sum_{i=1}^{\omega_d} B^{(d)}_i \Bigg \rfloor
$$

#A Simulation


#An Empirical Example from the GSS


#Discussion


#Conclusion


