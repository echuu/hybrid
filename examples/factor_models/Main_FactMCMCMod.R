

rm(list = ls())
#
# install.packages("tmvtnorm")  #used for the density of truncated multinormal
# install.packages("mvtnorm")
# install.packages("modeest")   # used for estimating mode
# install.packages("msm")       # used for generating univariate truncated normal
# install.packages("parallel")  # used for parallel computing
# install.packages("pscl")       # used for inverse gamma r.v.
# install.packages("MCMCpack")   # MCMCpack with MCMCfactanal
# install.packages("coda")
# install.packages("ars")


library(mvtnorm)
library(tmvtnorm)
library(modeest)
library(msm)
library(parallel)
library(pscl)
library(MCMCpack)
library(coda)
library(ars)

setwd("C:/Users/ericc/hybrid/examples")
source("GibbsFactMod.R")
source("ParAnalMod.R")
source("BetaMatrixBack.R")
source("LogPriorGivenKMod.R")
source("LogQsubK.R")
source("RJMCMCMod.R")
source("MyTrunNorm.R")
source("ParAnal.R")



# prior parameters for beta and Sigma

m = 4
k = 3
n = 15




# generate data when k = 3
set.seed(214)
# set.seed(213)
true.beta3 = matrix(0, nr = m, nc = 3)
true.beta3[, 1] =  c(rtnorm(n = 1, lower = 0, mean = -6), rnorm(n = m -1, mean = 9))
true.beta3[, 2] = c(0, rtnorm(n= 1, lower = 0, mean = 5), rnorm(n = m - 2, mean = 3))
true.beta3[, 3][3:m] = c(rtnorm(n= 1, lower = 0, mean = 2), rnorm(n = m - 3, mean = 5, sd = 18))
set.seed(123)
y3 = rmvnorm (n = n, mean= rep(0,m) , sigma = diag(m) + tcrossprod(true.beta3))
y3 = scale(y3, center = TRUE,  scale = TRUE)
y3.swap = y3[, sample(m)]

#
######################################################
######################################################


# initial beta values for parallel analysis
# only for beta's corresponding to k > 0

set.seed(1234)
k = 3
beta0 = matrix(0, nr = m, nc = k) # init value for beta0 used in Gibbs sampler
for(i in 1:k) {
    if(i == 1){
        beta0 [1, 1] = rtnorm(n = 1, mean = 0, sd =1, lower= 0)
    } else {
        lower.vec = c(rep(-Inf, i - 1), 0)
        beta0[i, 1:i] = rtmvnorm(n = 1, mean = rep(0, i), sigma = diag(i),
                                 lower = lower.vec)
    }

}
for(i in (k+1):m) {
    beta0[i,] = rmvnorm(n = 1, mean = rep(0, k), sigma = diag(k))
}

nu = 2.2
C0 = 1
s = sqrt(.1/2.2)

thin    = 20
burnin  = 1000
seed    = 123
No.Iter = 10000
Sigma0  = rep(1, m) ## init value for Sigma0 used in Gibbs sampler


obj = factGibbsMod(y = y3, k = k, beta0, rep(1, m), nu, C0, s,
                   No.Iter, burnin, thin, seed)

obj_swap = factGibbsMod(y = y3.swap, k = k, beta0, rep(1, m), nu, C0, s,
                   No.Iter, burnin, thin, seed)

str(obj)


### process the samples from factGibbsMod() function

formatBeta = function(x) {
  x[lower.tri(x, diag = TRUE)]
}

beta_samps = do.call(rbind, lapply(obj$beta, formatBeta))
omega_samps = do.call(rbind, obj$Sigma)
theta = cbind(beta_samps, omega_samps)









