#
# Example 2 in Chernozhukov-Hong (2003)
# Quantile IV Regression
# 
# David Kohns
# 13/03/22
#
rm(list=ls())
#
library(ggplot2); library(grid)
library(rstan)
library(rstanarm)
library(parallel)
library(doParallel)
options(mc.cores = parallel::detectCores())
numCores <- 4
registerDoParallel(numCores)

# Save Graphs:
saveG <- FALSE
#
set.seed(8981)
## ======== Generate Data ======== ##
K <- 100
N <- 500
L <- 100
beta <- matrix(rep(0,(K-1)))
beta[1:5,1]<- 2.5

# Now create a sample
X <- matrix(rnorm((K-1)*N),ncol=(K-1))
#X <- cbind(1,X)
Z <- cbind(1,scale(X))
epsilon <- rnorm(N)
#
sigmaD <- numeric(N)
for(i in 1:N){
  sigmaD[i] <- (sum(X[i,1:ncol(X)]))/4
}

# Heterskedastic errors
u <- sigmaD*epsilon
alpha <- 0
#
Y <- alpha + (X)%*%beta + matrix(u)
quantile(u,0.05)

## ========= Set Priors ========= ##
mu0 <- matrix(rep(0,K-1))
p0<- 15
scale_global<- p0/(K-p0)*1/sqrt(N)
scale_global <- 1 # Uncomment this for the normal rhs
nu_global <- 1
nu_local <- 1 # originally at 1
slab_scale <- 2 # originally at 2
slab_df <- 4 
lambda <- 1/(N^(1/4))

## ========= Compile Stan Model ========= ##
setwd("/Users/dk/Library/CloudStorage/OneDrive-Heriot-WattUniversity/Education/Non-chapter projects/Bayesian GMM")
writeLines(readLines("qiv_try_hs_naive.stan"))

test <- stan(file = 'qiv_try_rhs_gallant.stan',
             data = list(N=N,K=(K-1),L=K,mu=as.numeric(mu0),y=as.vector(scale((Y))),X=scale(X),Z=Z,tau=0.1,scale_global = scale_global,
                         nu_global = nu_global, nu_local = nu_local, slab_scale = slab_scale, slab_df = slab_df, a = 40, lambda_gallant = lambda),
             seed = 123456,
             control=list(adapt_delta=0.99, max_treedepth=15),
             warmup = 200,
             iter = 400)

test
test1 <- extract(test)
test1 <-colMeans( test1$beta)
test2 <- test1$beta[,6]
hist(test2)
AutoCorrelation <- acf(test2, plot = FALSE)
plot(AutoCorrelation, main = "Freeny's Quarterly Revenue Time Series ACF")
plot(test2)
hist(test2[10001:20000,4])
#END

library(R.matlab)
writeMat('/Users/dk/Library/CloudStorage/OneDrive-Heriot-WattUniversity/Education/Non-chapter projects/Bayesian GMM//k20_hmc.mat',scores_all=test2)
betatrue <- beta
betatrue[1] <- quantile(u,0.05) 
bias = sd(betatrue-test1)
writeMat('/Users/dk/Library/CloudStorage/OneDrive-Heriot-WattUniversity/Education/Non-chapter projects/Bayesian GMM//k20_hmc_bias.mat',scores_all=bias)