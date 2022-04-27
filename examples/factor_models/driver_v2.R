


set.seed(214)
m  = 4       # number of columns in y, number of rows in beta
k0 = 3       # number of columns in beta (number of factors)
n  = 500     # number of rows in y (number of observations)

tau = 0.25
omega0 = tau * diag(1, m) # identity matrix
beta0 = matrix(0, m, k0)
beta0[,1] = c(rtnorm(n = 1, lower = 0, mean = -1), rnorm(n = m - 1, mean = 1))
beta0[,2] = c(0, rtnorm(n = 1, lower = 0, mean = 1), rnorm(n = m - 2, mean = 1))
beta0[,3][3:m] = c(rtnorm(n= 1, lower = 0, mean = 1), rnorm(n = m - 3, mean = 1, sd = 1))
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
obj = factGibbsMod(y = y, k = k0, beta_init, Sigma_init, nu, C0, s,
                   No.Iter, burnin, thin)

## see if mean of betas is close to true beta
beta_list = form_beta(obj, length(obj$beta))
mean_beta(beta_list)
beta0
## see if covariance estimation is good
sigma_list = form_sigma(obj, length(obj$Sigma))
mean_sigma(sigma_list)
1/n * t(y) %*% y


### prepare object that stores model information to be passed into hyb-ep ######
k        = k0                ## number of factors in the fitted factor model
p        = m*(m+1)/2         ## number of elements in lower cholesky factor
K_mm     = matrixcalc::commutation.matrix(m, m)
ind_mat  = matrix(1:m^2, m, m)
chol_ind = ind_mat[lower.tri(ind_mat, diag = TRUE)]


# (1b) initialize parameter values for gibbs sampler
beta_init = matrix(0, nr = m, nc = k)
for(r in 1:k) {
  if(r == 1){
    beta_init[1, 1] = rtnorm(n = 1, mean = 0, sd = 1, lower = 0)
  } else {
    lower.vec = c(rep(-Inf, r - 1), 0)
    beta_init[r, 1:r] = rtmvnorm(n = 1, mean = rep(0, r), sigma = diag(r),
                                 lower = lower.vec)
  }
}
for(r in (k+1):m) {
  beta_init[r,] = rmvnorm(n = 1, mean = rep(0, k), sigma = diag(k))
}
Sigma_init = diag(omega0) ## init value for Sigma0 used in Gibbs sampler

## Gibbs sampler for the regular dataset: y
obj = factGibbsMod(y = y, k = k, beta_init, Sigma_init, nu, C0, s,
                   No.Iter, burnin, thin, seed)


#### draw samples from the prior of phi to find the sample mean/covariance
set.seed(123)
J = 1000 # number of prior samples to draw

# generate priors draws of sigma
sigma_prior = sample_sigma(J, m, k) ## this will change as we vary k (# factors)

# obtain prior draws of phi by computing the lower chol factor of sigma
phi_list = lapply(sigma_prior, getLowerChol) # phi * phi' = sigma

# estimate the mean, covariance of the approximating mvn distribution
prior = estimate_params(phi_list)

z0 = list(y = y, S = t(y) %*% y,
          n = n, nu = nu, C0 = C0, s = s, k = k, m = m, p = p,
          K_mm = K_mm, chol_ind = chol_ind,
          mu_phi = prior$mu,
          sigma_phi = prior$sigma)

#### obtain the posterior samples of sigma
#### using the posterior draws of (beta, sigma)

## run Gibbs sampler to obtain sampler object that contains posterior samples
## of (omega, beta). use these to compute sigma, and then compute lower cholesky
## and store the vectorized cholesky factors row-wise
phi_samps = phi_post(obj, n_post, p) # (n_post x p)

u_df = hybrid::preprocess(data.frame(phi_samps), p, z0)
u = phi_samps[1,]

## compute hybrid-ep estimator
library(dplyr)
hybrid::hybml(u_df, z0, psi = psi, grad = grad, hess = hess)$logz  # -4721.604

## compute laplace approximation at MAP
MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))[1]
ustar = u_df[MAP_LOC,1:p] %>% unname() %>% unlist()
0.5*(p)*log(2*pi) - 0.5*hybrid::log_det(hess(ustar, z0)) - psi(ustar, z0)

## compute modified vanilla hybrid
hybrid::hybml_map(u_df)$logz

## compute vanilla ep
hybridml::hybml_const(u_df)$logz



