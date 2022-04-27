
## log posterior of x = (beta, Omega)
logpost = function(x, z) { # x = (lower.tri(beta), diag(Omega))
  # x is (m+d)-dimensional
  # beta elements are d-dimensional, d = (m-k) * k + k * (k+1) / 2
  # omega elements are m-dimensional
  d = z$d
  m = z$m
  k = z$k
  n = z$n
  y = z$y

  xbeta = x[1:d]
  # xomega = x[(d+1):(m+d)]

  ## recreate beta matrix
  b = matrix(0, m, k)
  b[lower.tri(b, diag = TRUE)] = xbeta

  ## create Omega matrix
  # Omega = diag(xomega)
  Omega = z$Omega
  Sigma = Omega + b %*% t(b)
  Sigma_inv = solve(Sigma)

  q_term = 0
  for (i in 1:n) {
    q_term = q_term + c(t(y[i,]) %*% Sigma_inv %*% y[i,])
  }
  ## (1) compute loglikelihood -- this is the same value as ll_omega, ll_beta
  lp1 = -n*m/2 * log(2*pi) - n/2 * log(det(Sigma)) - 1/2 * q_term

  ## (2) compute log prior of beta
  lp2 = lp_beta(xbeta, z)

  ## (3) compute log prior of omega
  ## lp3 = lp_omega(xomega, z)

  return(lp1 + lp2)
} # end logpost() function


## gradient of log posterior of x = (beta, Omega)
grad_logpost = function(x, z) {
  d = z$d
  m = z$m
  k = z$k
  xbeta = x[1:d]
  # xomega = x[(d+1):(m+d)]

  ## recreate beta matrix
  beta = matrix(0, m, k)
  beta[lower.tri(beta, diag = TRUE)] = xbeta

  ## create Omega matrix
  # Omega = diag(xomega)
  Omega = z$Omega
  return(c(grad_logpost_beta(beta, Omega, z)))
  # return(c(grad_logpost_beta(beta, Omega, z),   ## grad of log-post wrt beta
  #          grad_logpost_omega(beta, Omega, z))) ## grad of log-post wrt Omega
} # end grad_logpost() function


## hessian of log posterior of x = (beta, Omega) -> (m + d) x (m + d) matrix
hess_logpost = function(x, z) {

  d = z$d
  m = z$m
  k = z$k
  y = z$y
  n = z$n
  xbeta = x[1:d]
  # xomega = x[(d+1):(m+d)]

  ## recreate beta matrix
  beta = matrix(0, m, k)
  beta[lower.tri(beta, diag = TRUE)] = xbeta
  ## create Omega matrix
  # Omega = diag(xomega)
  Omega = z$Omega

  # stack the hessian matrices along the diagonal to form block diagonal matrix
  # h = magic::adiag(hess_logpost_beta(beta, Omega, z),
  #                  hess_logpost_omega(beta, Omega, z))

  h = hess_logpost_beta(beta, Omega, z)

  # computed 2nd order mixed partials
  # h[z$beta_rows, z$omega_rows] = dbeta_domega(beta, Omega, z)
  # h[z$omega_rows, z$beta_rows] = t(h[z$beta_rows, z$omega_rows])

  return(h)
} # end hess_logpost() function


log_density = function(x, params) {
  -psi(x, params)
}
formatBeta = function(x) {
  x[lower.tri(x, diag = TRUE)]
}
psi = function(x, params) { # log likelihood + log prior for beta
  -logpost(x, params)
}
grad = function(x, params) {
  -grad_logpost(x, params)
}
hess = function(x, params) { # (d x d) since the parameter is only beta
  -hess_logpost(x, params)
}

preprocess = function(samps, d, params) {
  psi_u = apply(samps, 1, psi, params = params)
  # (1.2) name columns so that values can be extracted by partition.R
  u_df_names = character(d + 1)
  for (i in 1:d) {
    u_df_names[i] = paste("u", i, sep = '')
  }
  u_df_names[d + 1] = "psi_u"

  # populate u_df
  u_df = cbind(samps, psi_u) # J x (D + 1)
  names(u_df) = u_df_names
  return(u_df)
}

############### ----------------------------------------------------------------

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
n  = 300     # number of rows in y (number of observations)

# generate data when k = 3
# set.seed(214) # take this out when running sims
true.beta3 = matrix(0, nr = m, nc = k0)
true.beta3[, 1] =  c(rtnorm(n = 1, lower = 0, mean = -6), rnorm(n = m -1, mean = 9))
true.beta3[, 2] = c(0, rtnorm(n= 1, lower = 0, mean = 5), rnorm(n = m - 2, mean = 3))
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
K_mm = matrixcalc::commutation.matrix(m, m)z

ind_mat = matrix(1:(m*k), m, k)
lowerind = ind_mat[lower.tri(ind_mat, diag = TRUE)]
ind_mat = cbind((lowerind-1) %% m, floor((lowerind-1) / m) ) + 1
omega_index = (1:m - 1) * (m + 1) + 1
beta_rows = 1:d
omega_rows = (d+1):(d+m)

i_vec = 1:k
diag_index = 1 + (i_vec-1) * (m+1) - i_vec * (i_vec - 1) / 2


### initial values for Gibbs sampler
# set.seed(1234)
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

Sigma0  = rep(1, m) ## init value for Sigma0 used in Gibbs sampler

z0 = list(y = y, n = nrow(y), nu = nu, C0 = C0, s = s, k = k, m = m, d = d,
          ind_mat = ind_mat, K_mk = K_mk, K_mm = K_mm, beta_index = lowerind,
          omega_index = omega_index, beta_rows = beta_rows,
          omega_rows = omega_rows, diag_index = diag_index,
          Omega = Sigma0)


### start gibbs sampler --------------------------------------------------------

thin    = 20
burnin  = 3000
seed    = 123
No.Iter = 10000
## pass in initial values: (beta0, Sigma0)
obj = gibbsBeta(y = y, k = k, beta, Sigma0, nu, C0, s,
                No.Iter, burnin, thin, seed)

obj$beta %>% tail
beta

## compute hybrid approximation for this
# D = d + m
D = d # lower triangular + diag elements of beta are the only parameters now
beta_samps = do.call(rbind, lapply(obj$beta, formatBeta))
# omega_samps = do.call(rbind, obj$Sigma)
# theta = data.frame(cbind(beta_samps, omega_samps))
theta = data.frame(beta_samps)
u_df = preprocess(theta, D, z0)
u_df %>% head

hybridml::hybml_const(u_df)$logz


MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))[1]
ustar = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()

G = -hess(ustar, z0)
invG = solve(G)
thetaNew = ustar + 1 * invG %*% grad(ustar, z0)

psi(ustar, z0)
psi(thetaNew, z0)

# bridge_lb = rep(-Inf, D)
# bridge_ub = rep(Inf, D)
# bridge_approx(as.matrix(theta), log_density, z0, bridge_lb, bridge_ub)

ustar = hybrid::globalMode(u_df, z0)

hybrid::hybml(u_df, z0, psi = psi, grad = grad, hess = hess, u_0 = ustar)



hybrid::hybml_map(u_df)$logz
xbeta = beta0[lower.tri(beta0, diag = TRUE)]



# laplace-map estimator
0.5*(D)*log(2*pi) - 0.5*log_det(hess(ustar, z0)) - psi(ustar, z0)
# vanilla hybrid estimator using the global mode obj functions from hyb-ep
hybml_map(u_df)$logz
# vanilla hybrid estimator from aistat paper
hybridml::hybml_const(u_df)$logz







