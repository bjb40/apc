---
author:
	- name: Bryce Bartlett
	- name: Stephen Vaisey
title: Approximating Age-Period-Cohort Estimates by Averaging Multiple Models with Varying Window Restrictions
bibliography: citations/apc.bib
csl: citations/asa-mod.csl
---

#Abstract



Key words: Computational Sociology, Methods, Life Course, Bayesian, Culture

#Introduction


#Background

Both classic and modern studies have proposed different effects across three social dimensions of time: age, period, and cohort [@ryder_cohort_1965]. These dimensions of time operate on individuals in different ways as they flow through social space. Age effects are driven either by the common experiences associated with the biological process of aging [@jackson_biological_2003], or persistent age-structuring institutions, like high school or retirement [@mortimer_government_2003; @waite_constrained_2014].  Period effects are responses of everyone to contemporaneous social experiences, like recessions or wars [@lam_is_2014]. Finally, cohort effects are socialization effects. They are broad-based historical events that stick to the populations experiencing them, even after the event has long since ended [@vaisey_cultural_2016]. In contemporary research, cohorts are virtually always defined in terms of birth year.

These three dimensions of time sit at the center of some recurrent debates. For example, cultural sociologists are divided over whether the main driver of culture is a period process (cultural fragmentation) or a cohort process (acquired dispositions) [@vaisey_cultural_2016]. Clinical reasearchers and biodemographers are looking for biological age---indicators of biological deterioration---but argue about confounding from cohort changes [@jackson_biological_2003]. And yet others are concerned about separating long-term and short-term impacts of important events, like the Great Recession [@burgard_effects_2015].

Unfortunately, it is not trivial to estimate the unique effects of age, period, and cohort to find evidence that may settle these questions and debates. The fundamental problem is one of statistical identification: these three dimensions are linearly dependent. Two of the dimensions define the third. If a researcher knows an individual's age, and the year of the survey, cohort is also defined. Because these variables are exactly colinear, they are not estimable using classic statistical techniques. A number of solutions have been proposed to address this conundrum, these include traditional methods like block/window constraints [@glenn_strategies_2005], and new methods, including statistical transformation (the "intrinsic estimator"), and random effects models [@yang_age-period-cohort_2013]. 

The past several years have seen a resurgence in debates surrounding these issues [@luo_block_2016; @luo_assessing_2013; @bell_hierarchical_2017]. Perhaps the most common argument is that the "constraints", *i.e.* assumptions used to break the APC identity, produce unknown biases into the models. In reality, this presents an extreme case of multicollinearity, a well-known, though difficult problem in statistical theory [@wooldridge_introductory_2009, at p. 95-98], with numerous cautions in applied texts that separating highly colinear effects is difficult [@camm_essentials_2016, at p. 328]. From a Bayesian perspective, collinearity and multicollinearity significantly reduce the ability of the data to supply information to produce the estimates, potentially introducing sensitivity to the prior [@gelman_bayesian_2014, at p. 305-306], including assumtpions of linearity. To put it another way, the critique of APC models in general, and window constraints in particular [@luo_block_2016], is that different sets of assumptions lead to different (and inconsistent) results. But, what if Instead of focusing on finding a *best* fitting model, we use a simple nonlinear approach (block modelling), and hundreds of assumptions? BMA allows us to do just that: run hundreds of models with differeing assumptions, and then combine them, weighted by their fit with the data.

#Bayesian Model Averaging

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

where $\beta$ is an estimated effect, $\lambda$, $\rho$ and $\kappa$ index unique values for age, period, and cohort, and $A$, $P$, and $C$ stand for matricies of dummy variable series for age, period, and cohort. This model is unidentified, because, as with the continuous case, the dummy variables in any two of the matricies above fully condition the third matrix. In other words, the indicator variable in $C$ is a function of the indicators in $A$ and $P$ (in mathematical terms,the probability that any given cohort dummy variable is one or zero is exctly dependent on the values of A and P, so that $p(C_{\kappa}=0|A_{\lambda},P_{\kappa}$) is always exaclty zero, or exactly one, depending on the values of $A$ and $P$).

We can break this dependencey, however, by transofrming $A$, $P$, and $C$---preferably without resulting to some prespecified arbitrary set of constraints [@gelman_bayesian_2014, at p. 366]. How do the window constraints break the linear dependency? By way of example, we can construct an age dummy variable series where $A$ is sliced into two groups based on some cut-point, say $g$, so that ...

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


If the vector $G$ and $\lambda$ have the following properties: $max(G) = max(a)$, $min(G) < min(a)$, $\lambda \in  \{1,2, ... {\gamma-1}\}$, and $\gamma>2$, then $G$ can describe any posible sets of window restrictions for age. Generalizing cross all dimensions of APC, permuting three similar vectors (say $G^{(d)}$) will describe any model for any possible window constraints detailed in equation 1. Accordingly, all permutations of G, as defined above, constitute the model space of window constraints ($\mathscr{M}$). For any given set of data, $\mathscr{M}$ is finite, but it can be quite large. For example, 10 unique ages, periods, and cohorts present ```__``` possible models. ```describe how``` Of these```__``` models, only 1 is not estimable because of perfect colinearity. This target model is (theoretially) the least biased, although it is not the most parsimonious model. The question is how to best use information from some subset of possible models in $\mathscr{M}$ to estimate unbiased APC effects. Bayesian Model Averaging (BMA) provides a straightforward way to combine models.

The theoretical backdrop of BMA applies to this curcumstance quite well. In principal, there is no one true (or best) model; instead estimates are conditional on models from the modeling space ($\mathscr{M}$), and have  a posterior distribution, which is calculated as a weighted average of all models [@raftery_bayesian_1995, p. 144-145].  It has been applied in diverse areas from weather forecasting to biology to social science [@fragoso_bayesian_2015] .

##Bayesian Model Averaging (BMA) and the MC3 Algorithm 

Markov Chain methods for Bayesian Model Averaging (BMA) provides a sensible way to sample over a subset of continually better-fitting models and combine their estimates to produce an approximation of APC effects. The next two sections describe the MC3 agorithm developed for BMA, and outline the unique implementation of MC3 for this particular set of models, drawing from the Dirichelet distribution. 

##Using the Dirichlet Distribution to Sample Window Groups (G)

We simulate $G$ using two sets of nuisance parameters. For each dimension ($d$) of APC, the window breaks, $G^{(d)}$ of equaiton 2 are decomposed into (1) a cumulative sum from a Simplex for each dimension, with the same length of unique elements in $d$ ($B^{(d)}$), and (2) a scalar integer, $w^d$. $G^{(d)}$ is simply the product of $B$ times $w$. We sample $w$ from the uniform distibution, as follows:

$$
w^{(d)} \sim U(2,max(T^{(d)}))
$$

Where $\tau_{d}$ is the index number for the APC effects (the $\lambda$, $\rho$, or $\kappa$ of equaiton 1) $T$ is the continuous vector of values (the $A$,$P$, or $C$). We sample the break-points using a cumulative sum from a dirichelet distribution:

$$
B^{(d)} \sim Dir(\alpha_{\tau_d})
$$

The $\alpha$ for the Dirichelet distribution are for each of the vunique values in the dimension ($d$) for APC. The Dirchelet distribution draws a random projection on the standard simplex, *i.e.* a vector of weights $B$ weights between 0 and 1 which sum to 1. Multiplying the vector of weights ($B$) by a scalar $w$ provides a set of numbers which sum to the scalar value $w^{(d)}$. We transform $B$ and $w$ into matrix $G$ by selecting the unique floor rounded values of the product of the cumulative sums of $B$:

$$
G^{(d)}_b = \Bigg \lfloor w^{(d)} \sum_{i=1}^{\tau_d} B^{(d)}_i \Bigg \rfloor;
where G_b > G_{b-1} and b \subseteq \tau_d
$$

Where $G$ is a vector of the type described in equation 2, and $d$ indexes $G$ and is a subset of the unique values in $T$.

#A Simulation


#An Empirical Example from the GSS


#Discussion


#Conclusion


