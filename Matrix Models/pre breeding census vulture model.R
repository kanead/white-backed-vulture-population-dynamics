# -------------------------------------------------- 
# MATRIX POPULATION MODEL FOR AFRICAN WHITE BACKED VULTURES 
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
bp <- 0.85 # breeding propensity
cs <- 1 # clutch size
hs <- 0.76 # hatching success 
fs <- 0.6 # fledging success 
f1 <- bp * (cs/2) * hs * fs # divide by 2 to get females only

# survival 
s0 <- 0.42 # first year survival 
s1Kr <- 0.82 # juvenile survival Kruger
s2Kr <- 0.89 # subadult survival Kruger
s3Kr <- 1.0 # adult survival Kruger

#--------------------------------------------
# KRUGER POST-BREEDING CENSUS
#--------------------------------------------
# survival this year is multiplied by fecundity next year because in this model
# the birds have to survive the year before they become breeders i.e. from 4 years old to 
# breeding age at 5 years old - s2Kr*f1

# create the matrix for Kruger
MKrPost <- c(0,0,0,0,s2Kr*f1,s3Kr*f1,
         s0,0,0,0,0,0,
         0,s1Kr,0,0,0,0,
         0,0,s1Kr,0,0,0,
         0,0,0,s2Kr,0,0,
         0,0,0,0,s2Kr,s3Kr
)
MKrPost <- matrix ((MKrPost), ncol=6, byrow = TRUE)

# previous function is wrapped up into pop.projection
popKrugerPost<-eigen.analysis(MKrPost, zero=TRUE)
popKrugerPost$lambda1

#--------------------------------------------
# KRUGER PRE-BREEDING CENSUS
#--------------------------------------------
MKrpre <- c(0,0,0,0,s0*f1,
          s1Kr,0,0,0,0,
          0,s1Kr,0,0,0,
          0,0,s2Kr,0,0,
          0,0,0,s2Kr,s3Kr)
MKrpre <- matrix ((MKrpre), ncol=5, byrow = TRUE)

# previous function is wrapped up into pop.projection
popKrugerPre<-eigen.analysis(MKrpre, zero=TRUE)
popKrugerPre$lambda1

#--------------------------------------------
# KZN SURVIVAL RATES
#--------------------------------------------
s1Kz <-  0.8601882 # juvenile survival KZN
s2Kz <- 0.5134050 # subadult survival KZN
s3Kz <- 0.5672604 # adult survival KZN
#--------------------------------------------
# KZN POST-BREEDING CENSUS
#--------------------------------------------
# create the matrix for Kruger
MKZPost <- c(0,0,0,0,s2Kz*f1,s3Kz*f1,
             s0,0,0,0,0,0,
             0,s1Kz,0,0,0,0,
             0,0,s1Kz,0,0,0,
             0,0,0,s2Kz,0,0,
             0,0,0,0,s2Kz,s3Kz
)
MKZPost <- matrix ((MKZPost), ncol=6, byrow = TRUE)

# previous function is wrapped up into pop.projection
popKZPost<-eigen.analysis(MKZPost, zero=TRUE)
popKZPost$lambda1

#--------------------------------------------
# KZN PRE-BREEDING CENSUS
#--------------------------------------------
MKZpre <- c(0,0,0,0,s0*f1,
          s1Kz,0,0,0,0,
          0,s1Kz,0,0,0,
          0,0,s2Kz,0,0,
          0,0,0,s2Kz,s3Kz)
MKZpre <- matrix ((MKZpre), ncol=5, byrow = TRUE)

# previous function is wrapped up into pop.projection
popKZPre<-eigen.analysis(MKZpre, zero=TRUE)
popKZPre$lambda1
#--------------------------------------------
# MEGAMATRIX FOR METAPOPULATION STRUCTURE
#--------------------------------------------
# Effective migration rates (dispersal * stage-specific survival)
gb <- 0.1 # Kruger to KZN
bg <- 0.9 # KZN to Kruger

A <- matrix(c(
  0, 0, 0, 0, s0*f1, 0, 0, 0, 0, 0,
  s1Kr, 0, 0, 0, 0, s1Kz*bg, 0, 0, 0, 0,
  0, s1Kr*(1-gb), 0, 0, 0, 0, s1Kz*bg, 0, 0, 0,
  0, 0, s2Kr*(1-gb), 0, 0, 0, 0, s2Kz*bg, 0, 0,
  0, 0, 0, s2Kr*(1-gb), s3Kr*(1-gb), 0, 0, 0, s2Kz*bg, s3Kz*bg,
  0, 0, 0, 0, 0, 0, 0, 0, 0, s0*f1,
  s1Kr*gb, 0, 0, 0, 0, s1Kz, 0, 0, 0, 0,
  0, s1Kr*gb, 0, 0, 0, 0, s1Kz*(1-bg), 0, 0, 0,
  0, 0, s2Kr*gb, 0, 0, 0, 0, s2Kz*(1-bg), 0, 0,
  0, 0, 0, s2Kr*gb, s3Kr*gb, 0, 0, 0, s2Kz*(1-bg), s3Kz*(1-bg)), nrow = 10, byrow = TRUE)

A
lambda(A)
