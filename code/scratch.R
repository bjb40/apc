rm(list=ls())

source('config~.R')

dat.f = "H:/projects/gender_egal/output/private~/subdat.csv"

if(file.exists(dat.f)){
  dat = read.csv(dat.f)
  #save object for runing on the server
  save(dat,file=paste0(outdir,'~dat.RData')) #~ keeps it from pushing to git
} else{
  load(paste0(outdir,'~dat.RData'))
}

dat = dat %>%
  mutate(happy=3-happy,
         fechld = 4-fechld,
         fefam = fefam-1,
         fepol = fepol-1,
         fepresch = fepresch-1,
         egal = fechld + fepresch + fefam + fepol,
         birthyear = year - age,
         female=sex-1,
         race=factor(race,labels=c('White','Black','Other')))

#limit to complete cases // inlcude educaton only as covariate
dat = dat %>%
  dplyr::select(egal,age,year,birthyear,female,race,educ) %>%
  rename(a=age,
         c=birthyear,
         p=year)

#load fit object -- should be chains
load(paste0(outdir,'empirical_res2017-8-26.RData'))

use.wt = win$wt

#calculate opval
y.mean = mean(dat$egal)
ytilde = do.call(cbind,lapply(allmods,function(x) x$ytilde))
ytilde.m = colMeans(ytilde)
ytilde.mwt = sample(ytilde.m,1000,replace=TRUE,prob=use.wt)

#calculate equaltiy across each