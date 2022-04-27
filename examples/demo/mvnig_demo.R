
library(dplyr)

#### ----- define functions needed for hybrid, hybrid-ep estimators ----- ####

# psi(): negative log posterior
psi = function(u, params) {
  p = params$p
  # extract beta, sigmasq from posterior sample
  beta    = unname(unlist(u[1:p])) # (p x 1)
  sigmasq = unname(unlist(u[p+1])) # (1 x 1)

  # compute log posterior
  logpost = with(params,
                 sum(dnorm(y, X %*% beta, sqrt(sigmasq), log = T)) +
                   mvtnorm::dmvnorm(beta, mu_beta, sigmasq * V_beta, log = T) +
                   log(MCMCpack::dinvgamma(sigmasq, shape = a_0, scale = b_0)))
  return(-logpost)
} # end of psi() function

# For hybrid-ep estimator -- define gradient, hessian functions
grad = function(u, params) { pracma::grad(psi, x0 = u, params = params) }
hess = function(u, params) { pracma::hessian(psi, x0 = u, params = params) }

# sample_beta(): sample from [ beta | sigmasq, y ]
sample_beta = function(sigmasq, params) {
  mvtnorm::rmvnorm(1, mean = params$mu_n, sigma = sigmasq * params$V_n)
}

# sample_post(): sample from [ beta, sigmasq | y ]
sample_post = function(J, params) {
  sigmasq_post = with(params, MCMCpack::rinvgamma(J, shape = a_n, scale = b_n))
  beta_post = t(sapply(sigmasq_post, sample_beta, params = params))
  data.frame(beta_post, sigmasq_post)
}

# calculate true marginal likelihood
logML = function(params) {
  with(params,
       a_0 * log(b_0) + lgamma(a_n) + 0.5 * hybrid::log_det(V_n) -
         0.5 * hybrid::log_det(V_beta) - n / 2 * log(2 * pi) - lgamma(a_0) -
         a_n * log(b_n))
} # end of logML() function


set.seed(123)
d = 20                     # dimension of (beta, sigmasq)
n = 200                    # sample size
J = 100                    # number of posterior samples
p       = d - 1            # dimension of beta
mu_beta = rep(0, p)        # prior mean for beta
V_beta  = diag(1, p)       # prior scaled precision matrix for beta
a_0     = 1                # prior shape for sigmasq
b_0     = 0.5              # prior scale for sigmasq

## initalize true parameters
sigmasq = 4
beta    = sample(-10:10, p, replace = T)

## generate data using true parameters
X = matrix(rnorm(n * p), n, p)                # (n x p) design matrix
eps = rnorm(n, mean = 0, sd = sqrt(sigmasq))  # (n x 1) vector of residuals
y = X %*% beta + eps                          # (n x 1) response vector

### compute posterior parameters
V_beta_inv = solve(V_beta)
V_n_inv = t(X) %*% X + V_beta_inv

V_n  = solve(V_n_inv)                                # (p x p)
mu_n = V_n %*% (t(X) %*% y + V_beta_inv %*% mu_beta) # (p x 1)
a_n =  a_0 + n / 2
b_n =  c(b_0 + 0.5 * (t(y) %*% y +
                        t(mu_beta) %*% V_beta_inv %*% mu_beta -
                        t(mu_n) %*% V_n_inv %*% mu_n))

# store prior, posterior parameters
params = list(V_beta = V_beta, V_beta_inv = V_beta_inv, mu_beta = mu_beta,
              a_0 = a_0, b_0 = b_0,
              V_n = V_n, V_n_inv = V_n_inv, mu_n = mu_n,
              a_n = a_n, b_n = b_n,
              y = y, X = X, n = n, d = d, p = p)


# sample from posterior
samps = sample_post(J, params)
# evaluate the negative log posterior
u_df = hybrid::preprocess(samps, d, params)

# compute vanilla hybrid estimator
hybrid::hybml_const(u_df)$logz

# compute hybrid-ep estimator
hybrid::hybml(u_df, params, psi, grad, hess)$logz

# compute ground truth
logML(params)


