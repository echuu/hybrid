


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


set.seed(214)
m  = 3       # number of columns in y, number of rows in beta
k0 = 2       # number of columns in beta (number of factors)
n  = 500     # number of rows in y (number of observations)


tau = 0.25
omega0 = tau * diag(1, m) # identity matrix
beta0 = matrix(0, m, k0)
beta0[,1] = c(rtnorm(n = 1, lower = 0, mean = -1), rnorm(n = m - 1, mean = 1))
beta0[,2] = c(0, rtnorm(n = 1, lower = 0, mean = 1), rnorm(n = m - 2, mean = 1))
y = rmvnorm(n = n, mean = rep(0, m) , sigma = omega0 + tcrossprod(beta0))
# y3 = scale(y3, center = TRUE,  scale = TRUE)
swap_order = sample(m)
y_swap = y[, swap_order]

beta0 # true value of beta

# for now, we let the initial values be the true parameter values
Sigma_init = diag(omega0)
beta_init = beta0


nu = 2.2
C0 = 1
s = sqrt(.1/2.2)

thin    = 20
burnin  = 2000
seed    = 123
No.Iter = 1e4    # number of iterations to run Gibbs sampler
n_post = No.Iter / thin

## this gives us a net: No.Iter/thin posterior samples

## Gibbs sampler for the regular dataset: y
obj = factGibbsMod(y = y, k = k, beta_init, Sigma_init, nu, C0, s,
                   No.Iter, burnin, thin, seed)

# verify the parameter estimates are reasonaby close to true value of beta
beta_list = form_beta(obj, length(obj$beta))
mean_beta(beta_list)
beta0



### prepare object that stores model information to be passed into hyb-ep ######
k = k0
m = nrow(beta_init)
d = (m-k) * k + k * (k+1) / 2

K_mk = matrixcalc::commutation.matrix(m, k)
K_mm = matrixcalc::commutation.matrix(m, m)

ind_mat = matrix(1:(m*k), m, k)
lowerind = ind_mat[lower.tri(ind_mat, diag = TRUE)]
ind_mat = cbind((lowerind-1) %% m, floor((lowerind-1) / m) ) + 1
omega_index = (1:m - 1) * (m + 1) + 1
beta_rows = 1:d
omega_rows = (d+1):(d+m)

ind_mat = matrix(1:m^2, m, m)
chol_ind = ind_mat[lower.tri(ind_mat, diag = TRUE)]


i_vec = 1:k
diag_index = 1 + (i_vec-1) * (m+1) - i_vec * (i_vec - 1) / 2

z0 = list(y = y, n = nrow(y), nu = nu, C0 = C0, s = s, k = k, m = m, d = d,
          ind_mat = ind_mat, K_mk = K_mk, K_mm = K_mm, beta_index = lowerind,
          omega_index = omega_index, beta_rows = beta_rows,
          omega_rows = omega_rows, diag_index = diag_index,
          chol_ind = chol_ind,
          S = t(y) %*% y)



### compute hybrid approximation
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
psi(ustar, z0)

beta_star = matrix(0, m, k)
beta_star[lower.tri(beta_star, diag = TRUE)] = ustar[1:d]
omega_star = ustar[(d+1):D]

logpost(ustar, z0)

lp_beta(ustar[1:d], z0) + lp_omega(omega_star, z0)


hybrid::hybml(u_df, z0, psi = psi, grad = grad, hess = hess)



# laplace-map estimator
0.5*(D)*log(2*pi) - 0.5*hybrid::log_det(hess(ustar, z0)) - psi(ustar, z0)
# vanilla hybrid estimator using the global mode obj functions from hyb-ep
hybrid::hybml_map(u_df)$logz
# vanilla hybrid estimator from aistat paper
hybridml::hybml_const(u_df)$logz


#### laplace approximation using EP

H_star = hess(ustar, z0)
H_star_inv = solve(H_star)
lambda_star = grad(ustar, z0)
b_star = H_star %*% ustar - lambda_star
m_star = H_star_inv %*% (H_star %*% ustar - lambda_star)

lb = unname(apply(theta, 2, min))
ub = unname(apply(theta, 2, max))
mu = m_star
cov_mat = as.matrix(Matrix::forceSymmetric(H_star_inv))

gauss_int = hybrid::ep(mu, cov_mat, lb, ub)

D / 2 * log(2 * pi) - 0.5 * hybrid::log_det(H_star) -
  psi(ustar, z0) + sum(lambda_star * ustar) -
  0.5 * t(ustar) %*% H_star %*% ustar +
  0.5 * t(m_star) %*% H_star %*% m_star + gauss_int







############# helper/misc functions --------------------------------------------


form_sigma = function(obj, n) {
  sigma_list = vector(mode = "list", length = n)
  for (i in 1:n) {
    beta = matrix(unlist(obj$beta[i]), m, k)
    omega = diag(unlist(obj$Sigma[i]))
    sigma_list[[i]] = beta %*% t(beta) + omega
  }
  return(sigma_list)
}

form_beta = function(obj, n) {
  beta_list = vector(mode = "list", length = n)
  for (i in 1:n) {
    beta = matrix(unlist(obj$beta[i]), m, k)
    beta_list[[i]] = beta
  }
  return(beta_list)
}

form_beta_beta = function(obj, n) {
  beta_list = vector(mode = "list", length = n)
  for (i in 1:n) {
    beta = matrix(unlist(obj$beta[i]), m, k)
    beta_list[[i]] = beta %*% t(beta)
  }
  return(beta_list)
}

mean_sigma = function(sigma_list) {
  sigma_sum = matrix(0, m, m)
  for (i in 1:length(sigma_list)) {
    sigma_sum = sigma_sum + sigma_list[[i]]
  }
  return(sigma_sum/length(sigma_list))
}

mean_beta = function(beta_list) {
  beta_sum = matrix(0, m, k)
  for (i in 1:length(beta_list)) {
    beta_sum = beta_sum + beta_list[[i]]
  }
  return(beta_sum/length(beta_list))
}

mean_beta_beta = function(beta_list) {
  beta_sum = matrix(0, m, m)
  for (i in 1:length(beta_list)) {
    beta_sum = beta_sum + beta_list[[i]]
  }
  return(beta_sum/length(beta_list))
}


















