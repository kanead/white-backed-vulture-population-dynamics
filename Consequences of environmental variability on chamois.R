# Matrix Models for Population Management & Conservation
# 2014
# CNRS CEFE workshop, Montpellier, France
# Jean-Dominique Lebreton, Olivier Gimenez, Dave Koons

# EXERCISE 4: Consequences of environmental variability on chamois in the Bauges

rm(list=ls(all=TRUE))   # clears the R memory, which is sometimes useful

################################################################################
# Piece 1 of code: Deterministic model for a constant environment              #
################################################################################
library(MASS) # an R package needed for some matrix functions
library(popbio)

# First, define the demographic parameters
s0 <- 0.66
s1 <- 0.897
s210 <- 0.962
s11 <- 0.733
#s11 <- 0.962
f <- 1/2
m2 <- 0.66
m <- 0.92

# Next, create the 11x11 deterministic matrix model for chamois
A <- matrix(c(
  0, s0*m2*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,
  s1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, s210, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, s210, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, s210, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, s210, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, s210, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, s210, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, s210, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, s210, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, s210, s11),nrow=11, byrow=T)

# By now we have seen how to calculate the long-term population growth rate by
# using an eigen-analysis of the matrix model. Thus, we know what we are doing
# when we call the 'lambda' function from the 'popbio' package to do it for us.
lambda(A)

#To conduct a perturbation analysis of the model, we just need to put the string
#of demographic parameters and matrix model in 'lists', then call the vitalsens
#function from the popbio package.
chamois.vr <- list(s0=0.66,s1=0.897,s210=0.962,s11=0.733,f=1/2,m2=0.66,m=0.92)

chamois.A <- expression(
  0, s0*m2*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,
  s1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, s210, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, s210, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, s210, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, s210, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, s210, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, s210, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, s210, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, s210, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, s210, s11)

llsenselas <- vitalsens(chamois.A,chamois.vr)
llsenselas

################################################################################
# Piece 2 of code: Stochastic model for a variable environment                 #
################################################################################
# We have not yet examined stochastic matrix models in R, and there are many
# ways to model stochasticity. Thus, rather than using the 'popbio' package,
# lets develop some raw code to get a feeling for how to conduct a stochastic
# matrix model analysis.

# Specify the number of simulations, time steps for each simulation, and the
# pseudo-extinction threshold
sims <- 500
tspan <- 1000
threshold <- 100

# Define demographic parameters that do not vary over time
f <- 1/2
m2 <- 0.66
m <- 0.92

# Storage place for per time step growth rates for eventual calculation of
# the stochastic growth rate
gr <- matrix(0,sims,tspan-1)

# Storage for indicators on each simulation determining whether or not the
# population ever dropped below the pseudo-extinction threshold.
ext_ind <- matrix(0,sims,1)
extr_event <- matrix(0,sims,1)

for (j in 1:sims){
  # Define vector of initial abundance
  n <- c(100,100,100,100,100,100,100,100,100,100,100)
  nstore <- matrix(0,tspan,1) # temporary storage of time-specific abundance
  nstore[1] <- sum(n)
  for(t in 2:tspan){
    X <- rbinom(1, 1, 1/15)
    s0 <- 0.66 - (2/3*0.66*X)
    s1 <- 0.897 - (1/2*0.897*X)
    s210 <- 0.962 - (1/2*0.962*X)
    s11 <- 0.733 - (1/2*0.733*X)
    A <- matrix(c(
      0, s0*m2*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,s0*m*f,
      s1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, s210, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, s210, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, s210, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, s210, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, s210, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, s210, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, s210, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, s210, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, s210, s11),nrow=11, byrow=T)
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
