# --------------------------------------------------
# Vulture Metapopulation pseudo-extinction probability 
# --------------------------------------------------
# useful link http://www.mbr-pwrc.usgs.gov/workshops/uf2016/

# clean everything first
rm(list=ls())

# load required packages
library(popbio)
library(diagram)

#--------------------------------------------
# PARAMETERS 
#--------------------------------------------
# fecundity calculation, (Gauthier & Lebreton (2004) Population models for Greater Snow Geese)
bp <- 0.8 # breeding propensity
cs <- 1 # clutch size
hs <- 0.76 # hatching success 
fs <- 0.6 # fledging success 
f1 <- bp * (cs/2) * hs * fs # divide by 2 to get females only

# survival 
s0 <- 0.42 # first year survival # this value should probably be modified to account for 
# lower adult survival in KZN
s1Kr <- 0.82 # juvenile survival Kruger
s2Kr <- 0.89 # subadult survival Kruger
s3Kr <- 1.0 # adult survival Kruger

#--------------------------------------------
# KRUGER PRE-BREEDING CENSUS
#--------------------------------------------
MKrpre <- c(0,0,0,0,s0*f1,
            s1Kr,0,0,0,0,
            0,s1Kr,0,0,0,
            0,0,s2Kr,0,0,
            0,0,0,s2Kr,s3Kr)
MKrpre <- matrix ((MKrpre), ncol=5, byrow = TRUE)

lambda(MKrpre)

#--------------------------------------------
# KZN SURVIVAL RATES
#--------------------------------------------
s1Kz <-  0.86 # juvenile survival KZN
s2Kz <- 0.51 # subadult survival KZN
s3Kz <- 0.57 # adult survival KZN
#--------------------------------------------

#--------------------------------------------
# KZN PRE-BREEDING CENSUS
#--------------------------------------------
MKZpre <- c(0,0,0,0,s0*f1,
            s1Kz,0,0,0,0,
            0,s1Kz,0,0,0,
            0,0,s2Kz,0,0,
            0,0,0,s2Kz,s3Kz)
MKZpre <- matrix ((MKZpre), ncol=5, byrow = TRUE)

lambda(MKZpre)
#--------------------------------------------
# MODELLING EXTINCTION PROBABILITIES
#--------------------------------------------
# Specify the number of simulations, time steps for each simulation, and the
# pseudo-extinction threshold
sims <- 500
tspan <- 50
threshold <- 20

# Define demographic parameters that do not vary over time
f1 <- 0.1824
gb <- 0.01 # migration rates from Kruger to KZN
bg <- 0.01 # migration rates from KZN to Kruger
s0 <- 0.42
s1Kz <- 0.86
s2Kz <- 0.51
s3Kz <- 0.57
s1Kr <- 0.82
# Storage place for per time step growth rates for eventual calculation of
# the stochastic growth rate
gr <- matrix(0,sims,tspan-1)

# Storage for indicators on each simulation determining whether or not the
# population ever dropped below the pseudo-extinction threshold.
ext_ind <- matrix(0,sims,1)
extr_event <- matrix(0,sims,1)

for (j in 1:sims){
  # Define vector of initial abundance
  # juvenile < 2 years = 9%, immature, 3-5 years = 24%, adult > 5 years = 67%
  n <- c(54,54,144,144,804,100,100,100,100,100)
  nstore <- matrix(0,tspan,1) # temporary storage of time-specific abundance
  nstore[1] <- sum(n)
  for(t in 2:tspan){
    X <- rbinom(1, 1, 1/15)
    s2Kr <- 0.89 - (1/2*0.89*X)
    s3Kr <- 1 - (1/2*1*X)
    A <- matrix(c(
      0,  0,  0,  0,  s0*f1,  0,  0,  0,  0, 0,
      s1Kr*(1-gb), 0, 0, 0, 0, s1Kz*bg, 0, 0, 0, 0,
      0, s1Kr*(1-gb), 0, 0, 0, 0, s1Kz*bg, 0, 0, 0,
      0, 0, s2Kr*(1-gb), 0, 0, 0, 0, s2Kz*bg, 0, 0,
      0, 0, 0, s2Kr*(1-gb), s3Kr*(1-gb), 0, 0, 0, s2Kz*bg, s3Kz*bg,
      0, 0, 0, 0, 0, 0, 0, 0, 0, s0*f1,
      s1Kr*gb, 0, 0, 0, 0, s1Kz*(1-bg), 0, 0, 0, 0,
      0, s1Kr*gb, 0, 0, 0, 0, s1Kz*(1-bg), 0, 0, 0,
      0, 0, s2Kr*gb, 0, 0, 0, 0, s2Kz*(1-bg), 0, 0,
      0, 0, 0, s2Kr*gb, s3Kr*gb, 0, 0, 0, s2Kz*(1-bg), s3Kz*(1-bg)), nrow = 10, byrow = TRUE)
    nnew <- A%*%n
    nstore[t] <- sum(nnew)
    gr[j,t-1] <- sum(nnew)/sum(n)
    n <- nnew
  }
  if(min(nstore) < threshold) ext_ind[j] = 1
  if(X==1) extr_event[j]=1
}

ln_lambda_s <- mean(log(gr))      # the stochastic population growth rate
Lambda_s <- exp(ln_lambda_s)      # the stochastic population growth rate on the
nb_extr_event<-sum(extr_event)
# non-logged scale
ext_prob <- mean(ext_ind)         # the pseudo-extinction probability

ln_lambda_s
Lambda_s
ext_prob
nb_extr_event