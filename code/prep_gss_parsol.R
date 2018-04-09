#Load previously subset gss for gender egal
#
#Bryce Bartlett

rm(list=ls())

source('config~.R')
library(parallel)
library(lme4)
library(apcwin) #my library

#previously cleaned data for gender egalitarian project
#dat.f = "H:/projects/gender_egal/output/private~/subdat.csv"
dat.f = 'H:/projects/proposal/r_study/output/private~/cxgss.RData'


if(file.exists(dat.f)){
  load(dat.f)
  #save object for runing on the server
  save(cleandat,file=paste0(outdir,'~dat_parsol.RData')) #~ keeps it from pushing to git
} else{
  load(paste0(outdir,'~dat_parsol.RData'))
}

dat = cleandat %>%
  filter(!is.na(birthyear) & birthyear<=1982 & birthyear>1911 & 
           inballot==1 & year>=1994) %>%
  mutate(race=factor(race,labels=c('White','Black','Other')))

#limit to complete cases // inlcude educaton only as covariate
dat = dat %>%
  dplyr::select(parsol,age,year,birthyear,female,race,educ)

t=nrow(dat)

dat = dat[complete.cases(dat),]

cat('Proportion non-missing:',nrow(dat)/t)
rm(cleandat)


dat = dat %>% 
  rename(a=age,
         p=year,
         c=birthyear)

samps = apcsamp(dat,
                dv='parsol',
                cores=3,
                chains=3,
                method='ml',
                samples=250)

eff = draw_effs(samps,tol=0.01); rm(samps)

pp = plot(eff)

dat$ca = dat$a - mean(dat$a)

hapc = lmer(parsol~ca+I(ca^2) + factor(p) + factor(c) +
              (1|p)+ (1|c),
            data=dat)

##########


##########33
#combine them for an example...

save.image(file=paste0(outdir,'empirical_parsol.RData'))


