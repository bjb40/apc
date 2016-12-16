---
author: Bryce Bartlett
date: 12/15/2016
title: Combining estimates across models
tags: machine learning, bayes, model selection, methods, model averaging, meta-analysis
---

#Overview

I have a number of models with different breaks in the dependent variable. How do I combine estimats sensibly across models?

#Summary

These seem pretty straightforward ways to combine the estimates: meta-analysis, machine learning (random forests?), and Bayesian Model Averaging.

##Meta-analysis

I haven't done these meata-analyses, but have seen them applied in the context of gene expression and medicine. My approach would be similar to a Bayesian meta-analysis of different samples (as in BDA). But, there seems to be no reason not to do a fully Bayesian Model Averaging. There are detractors, but this is a fairly straightforward strategy.

##Machine Learning

Perhaps the most-used machine learning algorithms for this sort of analysis is an ensemble method known as Random Forests, which are an extension of classification and Regression Tree Algorithms [@breiman_classification_2001; @breiman_random_2001]. These can fairly easily be extended to a hierarchical bayesian model based on a dirichelet prior [@taddy_bayesian_2015]. The basic point of these is to partition the variables so that you have classifications (combinations) of variables, and do this over many iterations using a boot-strapped sub-sample which is cross-validated and "voted" by fit to the remainder (testing) sample.

There are large and robust iplementations in R (xgboost and Random Forests) and Python (sklearn and xgboost), and at least one user-defined functon in stat (https://ideas.repec.org/c/boc/bocode/s457932.html)[https://ideas.repec.org/c/boc/bocode/s457932.html]. This is good for me long term, but doesn't fit Scott's vision of the "window" constraint problem.

##Bayesian Model Averaging

In principal, there is no one true (or best) model; instead quantities estimated conditional on a model have a posterior distribution which is a weighted average of all models [@rafferty_1995, p. 144-145]. Raferty specifically identifies a regression parameter as a true value. The largest issue is one of computational intractibility --- it becomes impossible to estimate all models, so there are a number of data reduction strategies (such as Occam's window).

This is an area of active research for Rafferty, and he has a link to both research and packages to implement this: (http://www.stat.washington.edu/raftery/Research/bma.html)[http://www.stat.washington.edu/raftery/Research/bma.html]. Though Gelman and Lynch appear to have a relativley dim view of this as a strategy.

This fits Steve's strategy really well, and there are a finite set of models, and potential extensions from a tree-style model to reduce computational burdens. Notably, this is potentially able to combine with ensemble (machine-learning) methods [@].

#Conclusion

Should just try to do both in a Bayesian context!