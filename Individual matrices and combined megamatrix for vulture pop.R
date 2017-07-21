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
s0 <- 0.42 # first year survival # this value should probably be modified to account for 
# lower adult survival in KZN
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
s1Kz <-  0.86 # juvenile survival KZN
s2Kz <- 0.51 # subadult survival KZN
s3Kz <- 0.57 # adult survival KZN
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
gb <- 0.04 # Disperal from Kruger to KZN
bg <- 0.04 # Disperal from KZN to Kruger
# resight proportions from data 
# Kruger origin birds - 75 within 701 outside ~ 10% stay within the park
# KZN origin birds - 22 within 201 outside ~ 10% stay within KZN

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

A
lambda(A)

#--------------------------------------------
# SENSITIVITY AND ELASTICITY ANALYSIS
#--------------------------------------------
# Conduct a sensitivity and elasticity analysis for the lower-level vital rates
# (i.e, those that make up the matrix elements) using the popbio package.
# Just put the vital rates in a list, and write the matrix as an expression

vulture.vr <- list(s0=0.42,f1=0.1938,
                   s1Kr=0.82,s2Kr=0.89,s3Kr=1,
                   s1Kz=0.86,s2Kz=0.51,s3Kz=0.57,
                   gb=0.1,bg=0.1)

VultureA <- expression(
  0, 0, 0, 0, s0*f1, 0, 0, 0, 0, 0,
  s1Kr*(1-gb), 0, 0, 0, 0, s1Kz*bg, 0, 0, 0, 0,
  0, s1Kr*(1-gb), 0, 0, 0, 0, s1Kz*bg, 0, 0, 0,
  0, 0, s2Kr*(1-gb), 0, 0, 0, 0, s2Kz*bg, 0, 0,
  0, 0, 0, s2Kr*(1-gb), s3Kr*(1-gb), 0, 0, 0, s2Kz*bg, s3Kz*bg,
  0, 0, 0, 0, 0, 0, 0, 0, 0, s0*f1,
  s1Kr*gb, 0, 0, 0, 0, s1Kz*(1-bg), 0, 0, 0, 0,
  0, s1Kr*gb, 0, 0, 0, 0, s1Kz*(1-bg), 0, 0, 0,
  0, 0, s2Kr*gb, 0, 0, 0, 0, s2Kz*(1-bg), 0, 0,
  0, 0, 0, s2Kr*gb, s3Kr*gb, 0, 0, 0, s2Kz*(1-bg), s3Kz*(1-bg))

# then apply the following popbio function
llsenselas <- vitalsens(VultureA,vulture.vr)
llsenselas


#--------------------------------------------
# MEGAMATRIX FOR METAPOPULATION STRUCTURE SWITCHED
#--------------------------------------------
# Effective migration rates (dispersal * stage-specific survival)
gb <- 0.1 # Disperal from Kruger to KZN
bg <- 0.1 # Disperal from KZN to Kruger

ASwtich <- matrix(c(
  0,  0,  0,  0,  s0*f1,  0,  0,  0,  0, 0,
  s1Kz*(1-bg), 0, 0, 0, 0, s1Kr*gb, 0, 0, 0, 0,
  0, s1Kz*(1-bg), 0, 0, 0, 0, s1Kr*gb, 0, 0, 0,
  0, 0, s2Kz*(1-bg), 0, 0, 0, 0, s2Kr*gb, 0, 0,
  0, 0, 0, s2Kz*(1-bg), s3Kz*(1-bg), 0, 0, 0, s2Kr*gb, s3Kr*gb,
  0, 0, 0, 0, 0, 0, 0, 0, 0, s0*f1,
  s1Kz*bg, 0, 0, 0, 0, s1Kr*(1-gb), 0, 0, 0, 0,
  0, s1Kz*bg, 0, 0, 0, 0, s1Kr*(1-gb), 0, 0, 0,
  0, 0, s2Kz*bg, 0, 0, 0, 0, s2Kr*(1-gb), 0, 0,
  0, 0, 0, s2Kz*bg, s3Kz*bg, 0, 0, 0, s2Kr*(1-gb), s3Kr*(1-gb)), nrow = 10, byrow = TRUE)

ASwtich
lambda(ASwtich)




#--------------------------------------------
# MEGAMATRIX TEST FOR AGE-SPECIFIC EMIGRATION/IMMIGRATION
#--------------------------------------------
# Effective migration rates (dispersal * stage-specific survival)
gb0 <- 0.05
bg0 <- 0.05
gb <- 0.02 # Disperal from Kruger to KZN
bg <- 0.02 # Disperal from KZN to Kruger
gbA <- 0.05
bgA <- 0.05
# resight proportions from data 
# Kruger origin birds - 75 within 701 outside ~ 10% stay within the park
# KZN origin birds - 22 within 201 outside ~ 10% stay within KZN

Amig <- matrix(c(
  0,  0,  0,  0,  s0*(1-gb0)*f1,  0,  0,  0,  0, s0*bg0,
  s1Kr*(1-gb), 0, 0, 0, 0, s1Kz*bg, 0, 0, 0, 0,
  0, s1Kr*(1-gb), 0, 0, 0, 0, s1Kz*bg, 0, 0, 0,
  0, 0, s2Kr*(1-gb), 0, 0, 0, 0, s2Kz*bg, 0, 0,
  0, 0, 0, s2Kr*(1-gb), s3Kr*(1-gbA), 0, 0, 0, s2Kz*bg, s3Kz*bgA,
  0, 0, 0, 0, 0, 0, 0, 0, 0, s0*(1-bg0)*f1,
  s1Kr*gb, 0, 0, 0, s0*gb0, s1Kz*(1-bg), 0, 0, 0, 0,
  0, s1Kr*gb, 0, 0, 0, 0, s1Kz*(1-bg), 0, 0, 0,
  0, 0, s2Kr*gb, 0, 0, 0, 0, s2Kz*(1-bg), 0, 0,
  0, 0, 0, s2Kr*gb, s3Kr*gbA, 0, 0, 0, s2Kz*(1-bg), s3Kz*(1-bgA)), nrow = 10, byrow = TRUE)

Amig
lambda(Amig)

#--------------------------------------------
# MEGAMATRIX FOCUSSED ON KZN ONLY
#--------------------------------------------
# Effective migration rates (dispersal * stage-specific survival)
gb <- 0 # Disperal from Kruger to KZN
bg <- 0.01 # Disperal from KZN to Kruger

ASwtich <- matrix(c(
  0,  0,  0,  0,  s0*f1,  0,  0,  0,  0, 0,
  s1Kz*(1-bg), 0, 0, 0, 0, s1Kr*gb, 0, 0, 0, 0,
  0, s1Kz*(1-bg), 0, 0, 0, 0, s1Kr*gb, 0, 0, 0,
  0, 0, s2Kz*(1-bg), 0, 0, 0, 0, s2Kr*gb, 0, 0,
  0, 0, 0, s2Kz*(1-bg), s3Kz*(1-bg), 0, 0, 0, s2Kr*gb, s3Kr*gb,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  s1Kz*bg, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, s1Kz*bg, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, s2Kz*bg, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, s2Kz*bg, s3Kz*bg, 0, 0, 0, 0, 0), nrow = 10, byrow = TRUE)

ASwtich
lambda(ASwtich)

