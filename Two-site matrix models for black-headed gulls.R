# Matrix Models for Population Management & Conservation
# 2014
# CNRS CEFE workshop, Montpellier, France
# Jean-Dominique Lebreton, Olivier Gimenez, Dave Koons

# EXERCISE 3: Two-site matrix models for black-headed gulls

rm(list=ls(all=TRUE))   # clears the R memory, which is sometimes useful

################################################################################
# Piece 1 of code                                                              #
################################################################################
library(MASS) # an R package needed for some matrix functions
library(quadprog)
library(popbio)

# First, define the demographic parameters
# Good site
sg0 <- 0.4
sg1 <- 0.6
s <- 0.82  # adult survival common to both sites
fg <- 0.8  # half the mean clutch size (modeling females only)
pg2 <- 0.3 # proportion of adults that attempt to breed each year at age 2, etc.
pg3 <- 0.5
pg4 <- 0.7
pg5 <- 1
# Bad sites
sb0 <- 0.4
sb1 <- 0.5
fb <- 0.5  # half the mean clutch size at Bad sites
pb2 <- 0.5
pb3 <- 0.8
pb4 <- 1
pb5 <- 1
# Effective migration rates (juvenile dispersal * their survival)
s1gb <- 0.2 # Good to Bad
s1bg <- 0.3 # Bad to Good

# Next, create the two-site pre birth-pulse matrix model for black-headed gulls
A <- matrix(c(
  0, sg0*fg*pg2, sg0*fg*pg3, sg0*fg*pg4, sg0*fg*pg5, 0, 0, 0, 0, 0,
  sg1, 0, 0, 0, 0, s1bg, 0, 0, 0, 0,
  0, s, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, s, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, s, s, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, sb0*fb*pb2, sb0*fb*pb3, sb0*fb*pb4, sb0*fb*pb5,
  s1gb, 0, 0, 0, 0, sb1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, s, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, s, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, s, s), nrow = 10, byrow = TRUE)
A  
# Then use the following popbio function; quite easy!
lambda(A)

# The following code demonstrates the matrix algebra that
# is used 'behind the scences' of the lambda function in the popbio package.
rows <- dim(A)[1]
cols <- dim(A)[2]
eig <- eigen(A)           # eigenvalues of A
EigVecs <- eig$vectors    # eigenvectors of A
Lambdas <- Re(eig$values) # real number components of eigenvalues
Lambda <- max(Lambdas)    # long-term geometric rate of population growth
Lambda

pos <- which.max(Lambdas)  # finding the position of the dominant eigenvalue
w <- Re(eig$vectors[1:rows,pos]) # its associated right eigenvector
w
sad <- w/(sum(w))
sad <- round(sad,3) # scaled dominant right eigenvector: Stable Age Distribution
# In the following, the ginv function inverts a matrix by calculating the
# Moore-Penrose generalized inverse of a matrix.
V <- Conj(ginv(EigVecs))   # left eigenvector; NOTE this notation from H Caswell
v <- Re(t(t(V[pos,])))     # dominant left eigenvector
v
rv <- v/(sum(v))
rv <- round(rv,3)          # scaled to provide proportional Reproductive Values
sad
rv

# Conduct a sensitivity and elasticity analysis for the lower-level vital rates
# (i.e, those that make up the matrix elements) using the popbio package.
# Just put the vital rates in a list, and write the matrix as an expression
gull.vr <- list(sg0=0.4,sg1=0.6,s=0.82,fg=0.8,pg2=0.3,pg3=0.5,pg4=0.7,pg5=1,
                sb0=0.4,sb1=0.5,fb=0.5,pb2=0.5,pb3=0.8,pb4=1,pb5=1,s1gb=0.2,s1bg=0.3)

gull.A <- expression(
  0, sg0*fg*pg2, sg0*fg*pg3, sg0*fg*pg4, sg0*fg*pg5, 0, 0, 0, 0, 0,
  sg1, 0, 0, 0, 0, s1bg, 0, 0, 0, 0,
  0, s, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, s, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, s, s, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, sb0*fb*pb2, sb0*fb*pb3, sb0*fb*pb4, sb0*fb*pb5,
  s1gb, 0, 0, 0, 0, sb1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, s, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, s, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, s, s)

# then apply the following popbio function
llsenselas <- vitalsens(gull.A,gull.vr)
llsenselas


################################################################################
# Piece 2 of code                                                              #
################################################################################
ratio <- numeric()  # a storage bin for holding calculations of a ratio
# to be interpreted for question 5
n <- matrix(1,10,1) # a vector with an initial abundance of 1 individual per
# age and location
tspan <- 50
for (t in 1:tspan){
  n <- A%*%n        # %*% = matrix multiplication in R
  # note that we simply overwrite the abundance vector at each time step
  # since we have no need to store it for these questions. It is updating.
  Nbreedg <- pg2*n[2]+pg3*n[3]+pg4*n[4]+pg5*n[5]
  Nbreedb <- pb2*n[7]+pb3*n[8]+pb4*n[9]+pb5*n[10]
  ratio[t] <- Nbreedg/Nbreedb  # this we store
}
ratio

par(mar = c(5, 6, 4, 2))
plot(1:tspan,ratio,type="l",xlab=list("Time",cex=2),ylab=list("Ratio",cex=2))

################################################################################
# Piece 3 of code                                                              #
################################################################################
# Define adult dispersal probabilities
gb <- 0 # Good to Bad
bg <- 0 # Bad to Good

A <- matrix(c(
  0, sg0*fg*pg2, sg0*fg*pg3, sg0*fg*pg4, sg0*fg*pg5, 0, 0, 0, 0, 0,
  sg1, 0, 0, 0, 0, s1bg, 0, 0, 0, 0,
  0, s*(1-gb), 0, 0, 0, 0, s*bg, 0, 0, 0,
  0, 0, s*(1-gb), 0, 0, 0, 0, s*bg, 0, 0,
  0, 0, 0, s*(1-gb), s*(1-gb), 0, 0, 0, s*bg, s*bg,
  0, 0, 0, 0, 0, 0, sb0*fb*pb2, sb0*fb*pb3, sb0*fb*pb4, sb0*fb*pb5,
  s1gb, 0, 0, 0, 0, sb1, 0, 0, 0, 0,
  0, s*gb, 0, 0, 0, 0, s*(1-bg), 0, 0, 0,
  0, 0, s*gb, 0, 0, 0, 0, s*(1-bg), 0, 0,
  0, 0, 0, s*gb, s*gb, 0, 0, 0, s*(1-bg), s*(1-bg)), nrow = 10, byrow = TRUE)

ratio <- numeric()  # a storage bin for holding calculations of a ratio
# to be interpreted for question 5
n <- matrix(1,10,1) # a vector with an initial abundance of 1 individual per
# age and location
tspan <- 50
for (t in 1:tspan){
  n <- A%*%n        # %*% = matrix multiplication in R
  # note that we simply overwrite the abundance vector at each time step
  # since we have no need to store it for these questions. It is updating.
  Nbreedg <- pg2*n[2]+pg3*n[3]+pg4*n[4]+pg5*n[5]
  Nbreedb <- pb2*n[7]+pb3*n[8]+pb4*n[9]+pb5*n[10]
  ratio[t] <- Nbreedg/Nbreedb  # this we store
}

ratio