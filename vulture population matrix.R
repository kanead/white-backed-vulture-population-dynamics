# -------------------------------------------------- 
# Matrix Population Models for White-backed Vultures 
# --------------------------------------------------
# useful link http://www.mbr-pwrc.usgs.gov/workshops/uf2016/

# clean everything first
rm(list=ls())

# load required packages
library(popbio)
library(diagram)

# fecundity calculation, (Gauthier & Lebreton (2004) Population models for Greater Snow Geese)
bp <- 0.85 # breeding propensity
cs <- 1 # clutch size
hs <- 0.75 # hatching success 
fs <- 0.6 # fledging success 
fecundity <- bp * (cs/2) * hs * fs # divide by 2 to get females only
#--------------------------------------------
#               KRUGER
#--------------------------------------------
fsKr <- 0.42 # first year survival 
jsKr <- 0.7644702 # juvenile survival Kruger
ssKr <- 0.9256097 # subadult survival Kruger
asKr <- 0.9698704 # adult survival Kruger

# survival this year is multiplied by fecundity next year because in this model
# the birds have to survive the year before they become breeders i.e. from 4 years old to 
# breeding age at 5 years old 

# for Kruger

ssfKr <- ssKr * fecundity
asfKr <- asKr * fecundity

# create the matrix for Kruger
MKr <- c(0,0,0,0,ssfKr,asfKr,
         fsKr,0,0,0,0,0,
         0,jsKr,0,0,0,0,
         0,0,jsKr,0,0,0,
         0,0,0,ssKr,0,0,
         0,0,0,0,ssKr,asKr
)
MKr <- matrix ((MKr), ncol=6, byrow = TRUE)
colnames(MKr) <- c("babies","1yr olds","2yr olds","3yr olds","4yr olds","5yr olds")
MKr

# the process looks like the diagram with the nodes representing the age classes 
# Create population matrix
#
par(mfrow=c(1,1))

Numgenerations <- 6
DiffMat <- matrix(data = 0, nrow = Numgenerations, ncol = Numgenerations)
AA <- as.data.frame(DiffMat)
AA[[1,5]] <- "f[4]"
AA[[1,6]] <- "f[5]"
#
AA[[2,1]] <- "s[list(0,1)]"
AA[[3,2]] <- "s[list(1,2)]"
AA[[4,3]] <- "s[list(2,3)]"
AA[[5,4]] <- "s[list(3,4)]"
AA[[6,5]] <- "s[list(4,5)]"
AA[[6,6]] <- "s[list(5,5)]"

#
name <- c(expression(Age[0]), expression(Age[1]), expression(Age[2]),
          expression(Age[3]), expression(Age[4]), expression(Age[5]))
#
plotmat(A = AA, pos = 6, curve = 0.7, name = name, lwd = 2,
        arr.len = 0.6, arr.width = 0.25, my = -0.2,
        box.size = 0.05, arr.type = "triangle", dtext = 0.95,
        main = "Age-structured population model",
        relsize=0.97) 

# Population sizes at each age from Kruger

# Murn estimates 904 pairs of breeding adults 
# Assume an additional 0.3 immature and non-breeding birds per pair
# Remember, we're only modelling females 
additionalPopKr <- 904 * 0.3
totalPopKr <- 904*2 + additionalPopKr # = Murn's estimate
# divide additional population among the 5 non-adult categories 
# (Note: do we want some of these included among the adults as adult aged non-breeders?)
additionalPopKr / 5 / 2

nKR<-c(27, 27, 27, 27, 27, 904)
nKR<-matrix (nKR, ncol=1)
nKR

# previous function is wrapped up into pop.projection
popModelKr <- pop.projection(MKr,nKR,iterations=40)

# Calculate population growth rate and other demographic parameters from a projection matrix model
# using matrix algebra
eigen.analysis(MKr, zero=TRUE)

#--------------------------------------------
#               KZN
#--------------------------------------------
fsKZN <- 0.42 # first year survival 
jsKZN <-  0.9813814 # juvenile survival KZN
ssKZN <- 0.6928154 # subadult survival KZN
asKZN <- 0.5082008 # adult survival KZN

# survival this year is multiplied by fecundity next year because in this model
# the birds have to survive the year before they become breeders i.e. from 4 years old to 
# breeding age at 5 years old 

# for KZN

ssfKZN <- ssKZN * fecundity
asfKZN <- asKZN * fecundity

# create the matrix for KZN
MKZN <- c(0,0,0,0,ssfKZN,asfKZN,
          fsKZN,0,0,0,0,0,
          0,jsKZN,0,0,0,0,
          0,0,jsKZN,0,0,0,
          0,0,0,ssKZN,0,0,
          0,0,0,0,ssKZN,asKZN
)
MKZN <- matrix ((MKZN), ncol=6, byrow = TRUE)
colnames(MKZN) <- c("babies","1yr olds","2yr olds","3yr olds","4yr olds","5yr olds")
MKZN

# Population sizes at each age from KZN

# Rushworth estimates 319 pairs of breeding adults 
# Assume an additional 0.3 immature and non-breeding birds per pair
# Remember, we're only modelling females 
additionalPopKZN <- 319 * 0.3
totalPopKZN <- 319*2 + additionalPopKZN # < Rushworth's estimate
# divide additional population among the 5 non-adult categories 
# (Note: do we want some of these included among the adults as adult aged non-breeders?)
additionalPopKZN / 5 / 2
# this value is too low to reach the estimated 900 birds, instead we subtract the breeding population
# from the total pop estimate and divide the remainder up among the other 5 age categories
# Rushworth assumes there are between 800 and 900 birds in total in KZN, taking the 900 value
(900-319*2)/5/2


nKZN<-c(26, 26, 26, 26, 26, 319)
nKZN<-matrix (nKZN, ncol=1)
nKZN

# pop.projection function
popModelKZN <- pop.projection(MKZN,nKZN,iterations=40)

# Calculate population growth rate and other demographic parameters from a projection matrix model
# using matrix algebra
eigen.analysis(MKZN, zero=TRUE)

# Plot the data
# create panel plot to show the population trends of both populations side by side 
par(mfrow=c(1,2))
plot(popModelKr$pop.sizes, type="l", xlab = "year", ylab = "pop. size (females)", main = "Kruger")
plot(popModelKZN$pop.sizes, type="l", xlab = "year", ylab = "pop. size (females)", main = "KZN")
