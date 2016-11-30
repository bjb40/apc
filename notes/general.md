---
author: Bryce Bartlett
date: 11/29/2016
title: Block Model Reasoning
---

Two things to discuss -- (1) estimating the window constraints; and (2) combining the estimates.

# Window Constraints

The problem is not just an artifact of the data; it requires estimates of nonlinearities. There are three linearly related variables.

a. A similar problem is identifiable in terms of developmental changes. This is similar although easier to imagine because it is 2 dimensional instead of 4 dimensional. Assume that there are differential effects of age based on whether a person is an adult or not (adult defined as 18 and over). If you include two variables, age and age-18, the coefficients are linearly dependent and inestimable. There are well-defined ways to get out of this, however, including both block models and Piecewise linear growth models (e.g. Raudenbush 178-179). Though Raudenbush presents it in terms of random effects; it may also work in a fixed effects context.

b. Raudenbush also identifies the bias inherent in compositional changes (particularly for variances) when there are differences across age, and suggests a group-mean centering approach. This is very similar to Allison's "hybrid" random effects model; and *may* serve the same purpose of instrumenting out the time-varying nature of . The question here is how to mean center across AP, and C. The most obvious solution is some extension of the HAPC, and multiple different measures of age (as in a). 

c. One way to estimate the appropriate window constraints is to use an extension of Bayesian Change-point Models --- these identify nonlinearities as a hyper-parameter in a hierarchical model. Can be thought of sort of like a growth mixture model; except you don't identify the number of groups. It may be theoretically possible to shoe-horn the change-point models into a changeling framework, and then specify the number of windows ex ante, and use a model fit statistic to .

# Combining the estimates

These seem pretty straightforward ways to combine the estimates:

a. Some sort of meta-analysis. I haven't done one of these, but they are well-known (especially in medical literature). Scott also has a paper which combines regression results from a set of SEM equations. I believe this probably works out to be some sort of HAPC (based on recollection). Combining the window estimates is pretty straightforward. Specifically, I think you just transform the betas by multiplying with the window constraint, (including the associated covariance). (Scott and I have some simulations). Then all window constraints are on the same scale (*i.e.* one year). Another way to think of this is simply to construct the same design matrix to apply to all of the estimates (all relevant combinations with 1 unit increments across all APC categories). Then you just predict the outcomes, and the variances in the outcomes, and combine them using meta-analysis. Another place to look is genetics. I've seen some presentations, but they have well-known strategies to i.d. effects based on multiple testing. One issue here may be power.

b. Use the estimates to simulate and combine data data --- Scott and Trish do this in their recent paper. I haven't looked at the particulars, but this is similar to the Bayesian Model Averaging strategy.

c. Some Bayesian Model Averaging scheme. These are sensitive to prior specifications; however, because this isn't a variable selection model; it may be possible to come up with an appropriate prior specification. This is basically a weighted average of the estimates (which would have to be transformed, as above).

# Work flow

How do you want to do the logistics?

Start with Luo and Hodges's simulation, and make adjustments? (That would be the homerun, I suppose).