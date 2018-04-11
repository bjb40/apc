###
#this script gathers comparisons from tst_comparisons
source('config~.R')
#dependencies
library(apcwin) #custom package in dev
library(lme4) #random effects
library(merTools) #random effects display and sampling
library(parallel) #to speed up sampling
library(dplyr) #data manipulation
library(reshape2) # data manipulation
library(ggplot2) #plotting

load(paste0(outdir,'simdata_hapc/sim',lnum,'.RData'))