

library(mvtnorm)
library(tmvtnorm)
library(modeest)
library(msm)
library(parallel)
library(pscl)
library(MCMCpack)
library(coda)
library(ars)
library(dplyr)

setwd("C:/Users/ericc/hybrid/examples")
source("GibbsFactMod.R")
source("ParAnalMod.R")
source("BetaMatrixBack.R")
source("LogPriorGivenKMod.R")
source("LogQsubK.R")
source("RJMCMCMod.R")
source("MyTrunNorm.R")
source("ParAnal.R")

source("C:/Users/ericc/hybrid/examples/factor_model_density.R")
source("C:/Users/ericc/hybrid/examples/density_wrapper.R")

## generate true parameters, generate data using these parameters
m  = 3       # number of columns in y, number of rows in beta
k0 = 2       # number of columns in beta (number of factors)
n  = 20     # number of rows in y (number of observations)

# generate data when k = 3
# set.seed(214) # take this out when running sims
true.beta3 = matrix(0, nr = m, nc = k0)
true.beta3[, 1] =  c(rtnorm(n = 1, lower = 0, mean = -1), rnorm(n = m -1, mean = 1))
true.beta3[, 2] = c(0, rtnorm(n= 1, lower = 0, mean = 1), rnorm(n = m - 2, mean = 1))
# true.beta3[, 3][3:m] = c(rtnorm(n= 1, lower = 0, mean = 2), rnorm(n = m - 3, mean = 5, sd = 18))
set.seed(123)
y3 = rmvnorm (n = n, mean= rep(0,m) , sigma = diag(m) + tcrossprod(true.beta3))
# y3 = scale(y3, center = TRUE,  scale = TRUE)
y3.swap = y3[, sample(m)]


## fit model using a range for k = 2, 3, 4, 5
y = y3
beta = true.beta3
# Omega = diag(c(0.17, 0.05, 0.02, 0.02))
nu  = 2.2
C0  = 2
s   =  sqrt(.1/2.2)
C0  = 1; C_0 = 1;


## let kk be the number of factors we use in the fitted model

# k = ncol(beta)
k = 2
m = nrow(beta)
d = (m-k) * k + k * (k+1) / 2

K_mk = matrixcalc::commutation.matrix(m, k)
K_mm = matrixcalc::commutation.matrix(m, m)

ind_mat = matrix(1:(m*k), m, k)
lowerind = ind_mat[lower.tri(ind_mat, diag = TRUE)]
ind_mat = cbind((lowerind-1) %% m, floor((lowerind-1) / m) ) + 1
omega_index = (1:m - 1) * (m + 1) + 1
beta_rows = 1:d
omega_rows = (d+1):(d+m)

i_vec = 1:k
diag_index = 1 + (i_vec-1) * (m+1) - i_vec * (i_vec - 1) / 2

z0 = list(y = y, n = nrow(y), nu = nu, C0 = C0, s = s, k = k, m = m, d = d,
          ind_mat = ind_mat, K_mk = K_mk, K_mm = K_mm, beta_index = lowerind,
          omega_index = omega_index, beta_rows = beta_rows,
          omega_rows = omega_rows, diag_index = diag_index)


### initial values for Gibbs sampler
# set.seed(1234)
# beta0 = matrix(0, nr = m, nc = k) # init value for beta0 used in Gibbs sampler
# for(i in 1:k) {
#   if(i == 1){
#     beta0[1, 1] = rtnorm(n = 1, mean = 0, sd =1, lower= 0)
#   } else {
#     lower.vec = c(rep(-Inf, i - 1), 0)
#     beta0[i, 1:i] = rtmvnorm(n = 1, mean = rep(0, i), sigma = diag(i),
#                              lower = lower.vec)
#   }
#
# }
# for(i in (k+1):m) {
#   beta0[i,] = rmvnorm(n = 1, mean = rep(0, k), sigma = diag(k))
# }

Sigma0  = rep(1, m) ## init value for Sigma0 used in Gibbs sampler


#######

thin    = 20
burnin  = 25000
seed    = 123
No.Iter = 50000

thin    = 20
burnin  = 1000
seed    = 123
No.Iter = 10000

## pass in initial values: (beta0, Sigma0)
obj = factGibbsMod(y = y, k = k, beta, Sigma0, nu, C0, s,
                   No.Iter, burnin, thin, 123)
tail(obj$beta, 1)

obj_swap = factGibbsMod(y = y3.swap, k = k, beta, Sigma0, nu, C0, s,
                   No.Iter, burnin, thin, 123)
tail(obj_swap$beta, 1)

beta


beta_0 = matrix(unlist(obj$beta[1]), m, k)
omega_0 = diag(unlist(obj$Sigma[1]))
beta_0 %*% t(beta_0) + omega_0
diag(m) + tcrossprod(true.beta3)


## compute hybrid approximation for this
D = d + m
beta_samps = do.call(rbind, lapply(obj$beta, formatBeta))
omega_samps = do.call(rbind, obj$Sigma)
theta = data.frame(cbind(beta_samps, omega_samps))
u_df = preprocess(theta, D, z0)
u_df %>% tail

hybridml::hybml_const(u_df)$logz
hybrid::hybml_map(u_df)$logz

MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))[1]
ustar = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()

hybrid::hybml(u_df, z0, psi = psi, grad = grad, hess = hess)
hybrid::hybml(u_df, z0, psi = psi, grad = grad, hess = hess, u_0 = ustar)

xbeta = beta0[lower.tri(beta0, diag = TRUE)]



# laplace-map estimator
0.5*(D)*log(2*pi) - 0.5*hybrid::log_det(hess(ustar, z0)) - psi(ustar, z0)
# vanilla hybrid estimator using the global mode obj functions from hyb-ep
hybrid::hybml_map(u_df)$logz
# vanilla hybrid estimator from aistat paper
hybridml::hybml_const(u_df)$logz


