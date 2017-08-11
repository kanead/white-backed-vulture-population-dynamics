# -------------------------------------------------- 
# MATRIX POPULATION MODEL FOR AFRICAN WHITE BACKED VULTURES 
# --------------------------------------------------

# link http://www.mbr-pwrc.usgs.gov/workshops/uf2016/

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

# survival rate common to both 
s0 <- 0.42 # first year survival 
# this value should probably be modified to account for 
# lower adult survival in KZN

#--------------------------------------------
# KRUGER SURVIVAL RATES
#--------------------------------------------

s1Kr <- 0.82 # juvenile survival Kruger
s2Kr <- 0.89 # subadult survival Kruger
s3Kr <- 1.0 # adult survival Kruger

#--------------------------------------------
# KZN SURVIVAL RATES
#--------------------------------------------
s1Kz <-  0.86 # juvenile survival KZN
s2Kz <- 0.51 # subadult survival KZN
s3Kz <- 0.57 # adult survival KZN

#--------------------------------------------
# MEGAMATRIX TEST FOR AGE-SPECIFIC EMIGRATION/IMMIGRATION
#--------------------------------------------

# Effective migration rates 
gb0 <- 0.05 # 1st year migration Kruger to KZN
bg0 <- 0.05 # 1st year migration KZN to Kruger
gb <- 0.02 # 2nd year to 5th year migration Kruger to KZN
bg <- 0.02 # 2nd year to 5th year migration KZN to Kruger
gbA <- 0.05 # adult migration Kruger to KZN
bgA <- 0.05 # adult round




Amig <- matrix(c(
  0,  0,  0,  0,  s0*(1-gb0)*f1,  0,  0,  0,  0, s0*bg0,
  s1Kr*(1-gb), 0, 0, 0, 0, s1Kz*bg, 0, 0, 0, 0,
  0, s1Kr*(1-gb), 0, 0, 0, 0, s1Kz*bg, 0, 0, 0,
  0, 0, s2Kr*(1-gb), 0, 0, 0, 0, s2Kz*bg, 0, 0,
  0, 0, 0, s2Kr*(1-gb), s3Kr*(1-gbA), 0, 0, 0, s2Kz*bg, s3Kz*bgA,
  0, 0, 0, 0, s0*gb0, 0, 0, 0, 0, s0*(1-bg0)*f1,
  s1Kr*gb, 0, 0, 0, 0, s1Kz*(1-bg), 0, 0, 0, 0,
  0, s1Kr*gb, 0, 0, 0, 0, s1Kz*(1-bg), 0, 0, 0,
  0, 0, s2Kr*gb, 0, 0, 0, 0, s2Kz*(1-bg), 0, 0,
  0, 0, 0, s2Kr*gb, s3Kr*gbA, 0, 0, 0, s2Kz*(1-bg), s3Kz*(1-bgA)), nrow = 10, byrow = TRUE)

round(Amig,3)
lambda(Amig)


