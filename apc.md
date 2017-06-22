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

Unfortunately, it is not trivial to estimate the unique effects of age, period, and cohort to find evidence that may settle these debates. The fundamental problem is one of identification: these three dimensions are linearly dependent. Two of the dimensions define the third. 

#Description of the Method



##Windows: Piecewise Constant Functions



##The MC3 Algorithm and Bayesian Model Averaging

In principal, there is no one true (or best) model; instead quantities estimated conditional on a model have a posterior distribution which is a weighted average of all models [@rafferty_1995, p. 144-145]. Raferty specifically identifies a regression parameter as a true value. The largest issue is one of computational intractibility --- it becomes impossible to estimate all models, so there are a number of data reduction strategies (such as Occam's window).

##Dirichlet as a Window-Length Sampler


#A Simulation


#An Empirical Example from the GSS


#Discussion


#Conclusion


