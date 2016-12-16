---
author: Bryce Bartlett
date: 12/15/2016
title: Normal - DP Model
tags: bayes, methods, dirichlet process, stick-breaking
---

#Overview

This outlines a normal, DP hyperprior on the "window" constraints (such that individuals are assigned to different windows).

#Summary

##Likelihood

The normal model is fairly straightforward.

[[cut and paste from HLM]]

##Diricelet Process and Hyper-Priors

One alternative is to id a finite mixture hyperprior with the number of windows as an unknown parameter. The beauty of the DP process is it's direct applicability to the "windows" (and I especially like the stick breaking). Should probably try to implement both.

#Conclusion

Should just try to do both in a Bayesian context!