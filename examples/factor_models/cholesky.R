

Omega = omega0
sigma = omega0 + tcrossprod(beta0)
sigma_inv = solve(sigma)
phi = t(chol(sigma_inv)) # lower cholesky factor: sigma^{-1} = phi phi'
S = t(y) %*% y

beta_vec = beta0[lower.tri(beta0, diag = TRUE)]
ll_beta(beta_vec, z0)
x = c(beta_vec, diag(omega0))
x1 = unlist(unname(theta[1,]))

logpost(x, z0)
logpost(x1, z0)

ld_sigma = function(phi_vec) {
  phi = matrix(0, m, m)
  phi[lower.tri(phi, diag = TRUE)] = phi_vec
  sigma = solve(phi %*% t(phi))
  log(det(sigma))
}

grad_ld_sigma = function(phi_vec) {
  phi = matrix(0, m, m)
  phi[lower.tri(phi, diag = TRUE)] = phi_vec
  sigma = solve(phi %*% t(phi))
  dphi = - 2 * sigma %*% phi
  dphi[lower.tri(dphi, diag = TRUE)]
}

hess_ld_sigma = function(phi_vec, z) {
  m = z$m
  lower_ind = z$chol_ind
  phi = matrix(0, m, m)
  phi[lower.tri(phi, diag = TRUE)] = phi_vec
  sigma = solve(phi %*% t(phi))
  I_m = diag(1, m)
  K_mm = z$K_mm

  h = - 2 * (I_m %x% sigma) +
    (((t(phi) %*% sigma) %x% (sigma %*% phi)) %*% K_mm +
       t(K_mm) %*% ((sigma %*% phi) %x% (t(phi) %*% sigma))) +
    2 * ((t(phi) %*% sigma %*% phi) %x% sigma)
  h[lower_ind, lower_ind]
}


etr_sigma = function(phi_vec, z0) {
  phi = matrix(0, m, m)
  phi[lower.tri(phi, diag = TRUE)] = phi_vec
  matrixcalc::matrix.trace(phi %*% t(phi) %*% S)
}

grad_etr_sigma = function(phi_vec) {
  phi = matrix(0, m, m)
  phi[lower.tri(phi, diag = TRUE)] = phi_vec
  dphi = 2 * S %*% phi
  dphi[lower.tri(dphi, diag = TRUE)]
}

hess_etr_sigma = function(phi_vec, z) {
  m = z$m
  S = z$S
  lower_ind = z$chol_ind
  I_m = diag(1, m)
  h = 2 * I_m %x% S
  # return(h)
  h[lower_ind, lower_ind]
}


# ------------------------------------------------------------------------------

loglik = function(phi_vec, z) {
  # x is (m+d)-dimensional
  # beta elements are d-dimensional, d = (m-k) * k + k * (k+1) / 2
  # omega elements are m-dimensional
  d = z$d
  m = z$m
  k = z$k
  n = z$n
  S = z$S

  phi = matrix(0, m, m)
  phi[lower.tri(phi, diag = TRUE)] = phi_vec

  sigma_inv = phi %*% t(phi)

  ## (1) compute loglikelihood -- this is the same value as ll_omega, ll_beta
  lp1 = -n*m/2 * log(2*pi) + n/2 * log(det(sigma_inv)) -
    1/2 * matrixcalc::matrix.trace(sigma_inv %*% S)

  return(lp1)
}

grad_loglik = function(phi_vec, z) {
  m = z$m
  S = z$S
  n = z$n

  phi = matrix(0, m, m)
  phi[lower.tri(phi, diag = TRUE)] = phi_vec

  sigma_inv = phi %*% t(phi)
  sigma = solve(sigma_inv)

  ## compute gradient of log-determinant term
  dphi = - 2 * sigma %*% phi
  g1 = dphi[lower.tri(dphi, diag = TRUE)]

  ## compute gradient of trace term
  dphi = 2 * S %*% phi
  g2 = dphi[lower.tri(dphi, diag = TRUE)]

  - n/2 * g1 - 1/2 * g2
}

hess_loglik = function(phi_vec, z) {
  n = z$n
  m = z$m
  S = z$S
  phi = matrix(0, m, m)
  phi[lower.tri(phi, diag = TRUE)] = phi_vec
  sigma = solve(phi %*% t(phi))
  I_m = diag(1, m)
  K_mm = z$K_mm
  lower_ind = z$chol_ind

  h1 = - 2 * (I_m %x% sigma) +
    (((t(phi) %*% sigma) %x% (sigma %*% phi)) %*% K_mm +
       t(K_mm) %*% ((sigma %*% phi) %x% (t(phi) %*% sigma))) +
    2 * ((t(phi) %*% sigma %*% phi) %x% sigma)

  h2 = 2 * I_m %x% S

  h = -n/2 * h1 - 1/2 * h2
  h[lower_ind, lower_ind]
}


loglik(phi_vec, z0)
pracma::grad(loglik, x0 = phi_vec, z = z0)
grad_loglik(phi_vec, z0)

pracma::hessian(loglik, x0 = phi_vec, z = z0)
hess_loglik(phi_vec, z0)

pracma::hessian(ld_sigma, x0 = phi_vec)
hess_ld_sigma(phi_vec, z0)

pracma::hessian(etr_sigma, phi_vec, z0 = z0)
hess_etr_sigma(phi_vec, z0)

# ------------------------------------------------------------------------------

## log prior functions

lp_phi = function(x, z) {
  p = z$p
  sigma_phi = z$sigma_phi
  mu_phi = z$mu_phi
  c(-p/2 * log(2*pi) - 1/2 * log(det(sigma_phi)) -
    1/2 * t(x-mu_phi) %*% solve(sigma_phi) %*% (x-mu_phi))
}

grad_lp_phi = function(x, z) {
  c(solve(z$sigma_phi) %*% (z$mu_phi - x))
}


hess_lp_phi = function(x, z) {
  return(-solve(z$sigma_phi))
}

# log prior evaluation of lower cholesky factor
lp_phi(phi_vec, z0)
dmvnorm(phi_vec, z0$mu_phi, z0$sigma_phi, log = TRUE)

# gradient of log prior
grad_lp_phi(phi_vec, z0)
pracma::grad(lp_phi, x0 = phi_vec, z = z0)

# hessian of log prior
hess_lp_phi(phi_vec, z0)
pracma::hessian(lp_phi, x0 = phi_vec, z = z0)

all.equal(hess_lp_phi(phi_vec, z0),
          pracma::hessian(lp_phi, x0 = phi_vec, z = z0))

# ------------------------------------------------------------------------------

#### tmp definition of log like for different parametrization of phi

# loglik = function(phi_vec, z) {
#   # x is (m+d)-dimensional
#   # beta elements are d-dimensional, d = (m-k) * k + k * (k+1) / 2
#   # omega elements are m-dimensional
#   d = z$d
#   m = z$m
#   k = z$k
#   n = z$n
#   S = z$S
#
#   phi = matrix(0, m, m)
#   phi[lower.tri(phi, diag = TRUE)] = phi_vec
#
#   sigma = phi %*% t(phi)
#   sigma_inv = solve(sigma)
#
#   ## (1) compute loglikelihood -- this is the same value as ll_omega, ll_beta
#   lp1 = -n*m/2 * log(2*pi) + n/2 * log(det(sigma_inv)) -
#     1/2 * matrixcalc::matrix.trace(sigma_inv %*% S)
#
#   return(lp1)
# }
#
# phi_post = function(obj, J, p) {
#   samps = matrix(0, J, p) # store vectorized form of phi in matrix row-wise
#   for (j in 1:J) {
#     omega = diag(obj$Sigma[[j]])
#     beta = obj$beta[[j]]
#     sigma = omega + beta %*% t(beta)
#     # sigma_inv = solve(sigma)
#     phi = getLowerChol(sigma)
#     samps[j,] = phi[lower.tri(phi, diag = TRUE)]
#   }
#   return(samps)
# }
#
# sample_sigma = function(J, m, k) {
#   omega_samps = matrix(
#     pscl::rigamma(m * J, alpha = nu / 2, beta = nu * s^2 / 2), J, m
#   )
#
#   sigma_list = vector(mode = "list", length = J)
#   for (j in 1:J) {
#     omega = diag(omega_samps[j,])
#     beta = sample_beta(m, k)
#     sigma = omega + beta %*% t(beta)
#     sigma_list[[j]] = sigma ## store sigma for alternative param
#   }
#   return(sigma_list)
# }
#
# psi = function(x, params) {
#   -logpost_phi(x, params)
# }
# grad = function(x, params) {
#   # -grad_logpost_phi(x, params)
#   pracma::grad(psi, x, params = params)
# }
# hess = function(x, params) {
#   # -hess_logpost_phi(x, params)
#   pracma::hessian(psi, x, params = params)
# }


#### end of alternative parametrization of phi




## log posterior functions

logpost_phi = function(x, z) {
  loglik(x, z) + lp_phi(x, z)
}
grad_logpost_phi = function(x, z) {
  grad_loglik(x, z) + grad_lp_phi(x, z)
}
hess_logpost_phi = function(x, z) {
  hess_loglik(x, z) + hess_lp_phi(x, z)
}

## verify equality of numerical, analyic gradient, hessian
pracma::grad(logpost_phi, x0 = phi_vec, z = z0)
grad_logpost_phi(phi_vec, z0)

pracma::hessian(logpost_phi, x0 = phi_vec, z = z0)
hess_logpost_phi(phi_vec, z0)

## general wrapper functions that vanilla hybrid + hybrid-ep will both use
psi = function(x, params) {
  -logpost_phi(x, params)
}
grad = function(x, params) {
  -grad_logpost_phi(x, params)
}
hess = function(x, params) {
  -hess_logpost_phi(x, params)
}




