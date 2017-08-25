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
eigen.analysis(MKrpre)
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
eigen.analysis(MKZpre)

#--------------------------------------------
# METAPOPULATION MATRIX
#--------------------------------------------
f1 <- 0.1824 * 0.42 # fecundity * first year survival s0
gb0 <- 0.05 # migration rates from Kruger to KZN for ages < 1
bg0 <- 0.05 # migration rates from KZN to Kruger for ages < 1
gb <- 0.02 # migration rates from Kruger to KZN for ages > 1  
bg <- 0.02 # migration rates from KZN to Kruger for ages > 1
s1Kz <- 0.86 # juvenile survival KZN
s2Kz <- 0.51 # Subadult survival KZN
s3Kz <- 0.57 # adult survival KZN
s1Kr <- 0.82 # juvenile survival Kruger
s2Kr <- 0.89 # Kruger subadult survival
s3Kr <- 1 # Kruger adult survival

metaA <- matrix(c( 
  # Kruger                                            # KZN
  0,  0,  0,  0,  f1*(1-gb0),                          0,  0,  0,  0, f1*bg0,
  s1Kr*(1-gb), 0, 0, 0, 0,                            s1Kz*bg, 0, 0, 0, 0,
  0, s1Kr*(1-gb), 0, 0, 0,                            0, s1Kz*bg, 0, 0, 0,
  0, 0, s2Kr*(1-gb), 0, 0,                            0, 0, s2Kz*bg, 0, 0,
  0, 0, 0, s2Kr*(1-gb), s3Kr*(1-gb),                  0, 0, 0, s2Kz*bg, s3Kz*bg,
  
  0, 0, 0, 0, f1*gb0,                                  0, 0, 0, 0, f1*(1-bg0),
  s1Kr*gb, 0, 0, 0, 0,                                s1Kz*(1-bg), 0, 0, 0, 0,
  0, s1Kr*gb, 0, 0, 0,                                0, s1Kz*(1-bg), 0, 0, 0,
  0, 0, s2Kr*gb, 0, 0,                                0, 0, s2Kz*(1-bg), 0, 0,
  0, 0, 0, s2Kr*gb, s3Kr*gb,                          0, 0, 0, s2Kz*(1-bg), s3Kz*(1-bg)), 
  
  nrow = 10, byrow = TRUE)

lambda(metaA)



#--------------------------------------------
# SENSITIVITY AND ELASTICITY ANALYSIS
#--------------------------------------------
# Conduct a sensitivity and elasticity analysis for the lower-level vital rates
# (i.e, those that make up the matrix elements) using the popbio package.
# Just put the vital rates in a list, and write the matrix as an expression

vulture.vr <- list(f1 = 0.076608, # fecundity * first year survival s0 = 0.1824 * 0.42
                   gb0 = 0.05, # migration rates from Kruger to KZN for ages < 1
                   bg0 = 0.05, # migration rates from KZN to Kruger for ages < 1
                   gb = 0.02, # migration rates from Kruger to KZN for ages > 1  
                   bg = 0.02, # migration rates from KZN to Kruger for ages > 1
                   s1Kz = 0.86, # juvenile survival KZN
                   s2Kz = 0.51, # Subadult survival KZN
                   s3Kz = 0.57, # adult survival KZN
                   s1Kr = 0.82, # juvenile survival Kruger
                   s2Kr = 0.89, # Kruger subadult survival
                   s3Kr = 1) # Kruger adult survival

                   
sensA <- expression( 
                     # Kruger                                            # KZN
                     0,  0,  0,  0,  f1*(1-gb0),                          0,  0,  0,  0, f1*bg0,
                     s1Kr*(1-gb), 0, 0, 0, 0,                            s1Kz*bg, 0, 0, 0, 0,
                     0, s1Kr*(1-gb), 0, 0, 0,                            0, s1Kz*bg, 0, 0, 0,
                     0, 0, s2Kr*(1-gb), 0, 0,                            0, 0, s2Kz*bg, 0, 0,
                     0, 0, 0, s2Kr*(1-gb), s3Kr*(1-gb),                  0, 0, 0, s2Kz*bg, s3Kz*bg,
                     
                     0, 0, 0, 0, f1*gb0,                                  0, 0, 0, 0, f1*(1-bg0),
                     s1Kr*gb, 0, 0, 0, 0,                                s1Kz*(1-bg), 0, 0, 0, 0,
                     0, s1Kr*gb, 0, 0, 0,                                0, s1Kz*(1-bg), 0, 0, 0,
                     0, 0, s2Kr*gb, 0, 0,                                0, 0, s2Kz*(1-bg), 0, 0,
                     0, 0, 0, s2Kr*gb, s3Kr*gb,                          0, 0, 0, s2Kz*(1-bg), s3Kz*(1-bg)
                     )

# then apply the following popbio function
llsenselas <- vitalsens(sensA,vulture.vr)
llsenselas


#--------------------------------------------
# MODELLING EXTINCTION PROBABILITIES
#--------------------------------------------
# Specify the number of simulations, time steps for each simulation, and the
# pseudo-extinction threshold
sims <- 500 
tspan <- 100
threshold <- 20

# Define demographic parameters that do not vary over time
f1 <- 0.1824 * 0.42 # fecundity * first year survival s0
gb0 <- 0.05 # migration rates from Kruger to KZN for ages < 1
bg0 <- 0.05 # migration rates from KZN to Kruger for ages < 1
gb <- 0.02 # migration rates from Kruger to KZN for ages > 1  
bg <- 0.02 # migration rates from KZN to Kruger for ages > 1
s1Kz <- 0.86 # juvenile survival KZN
s2Kz <- 0.51 # Subadult survival KZN
s3Kz <- 0.57 # adult survival KZN
s1Kr <- 0.82 # juvenile survival Kruger

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
  # Kruger 1200 females 
  # < 2 = 108
  # 3-5 = 288 / 3 = 96
  # 5+ = 804
  # KZN 425 females 
  # < 2 = 38
  # 3-5 = 102 / 3 = 34
  # 5+ = 285
  
  n <- c(108,96,96,96,804,38,34,34,34,285)
  nstore <- matrix(0,tspan,1) # temporary storage of time-specific abundance
  nstore[1] <- sum(n)
  for(t in 2:tspan){
    X <- rbinom(1, 1, 1/15) # probabilty that poisoning occurs, here it's once ever 10 years
    s2Kr <- 0.89 - (0.427*0.89*X) # Kruger subadult survival equals KZN under poisoning
    s3Kr <- 1 - (0.43*1*X) # Kruger adult survival equals KZN under poisoning
    
# The following megamatrix matrix - A - has the Kruger matrix in top left 
# and KZN matrix in bottom right
# the diagonals are include the probability of emigration/immigration between sites
# that's why it's called a megamatrix
# take a look at the simpler matrices above constructed separately for each site
# to see where these values come from
    
    A <- matrix(c( 
      # Kruger                                            # KZN
      0,  0,  0,  0,  f1*(1-gb0),                          0,  0,  0,  0, f1*bg0,
      s1Kr*(1-gb), 0, 0, 0, 0,                            s1Kz*bg, 0, 0, 0, 0,
      0, s1Kr*(1-gb), 0, 0, 0,                            0, s1Kz*bg, 0, 0, 0,
      0, 0, s2Kr*(1-gb), 0, 0,                            0, 0, s2Kz*bg, 0, 0,
      0, 0, 0, s2Kr*(1-gb), s3Kr*(1-gb),                  0, 0, 0, s2Kz*bg, s3Kz*bg,
      
      0, 0, 0, 0, f1*gb0,                                  0, 0, 0, 0, f1*(1-bg0),
      s1Kr*gb, 0, 0, 0, 0,                                s1Kz*(1-bg), 0, 0, 0, 0,
      0, s1Kr*gb, 0, 0, 0,                                0, s1Kz*(1-bg), 0, 0, 0,
      0, 0, s2Kr*gb, 0, 0,                                0, 0, s2Kz*(1-bg), 0, 0,
      0, 0, 0, s2Kr*gb, s3Kr*gb,                          0, 0, 0, s2Kz*(1-bg), s3Kz*(1-bg)), 
      
      nrow = 10, byrow = TRUE)
    
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
ext_prob # the overall extinction probability 
nb_extr_event # number of extinction events over the model run