


## using the the prior samples of phi, estimate the mean and covariance of the
## approximating MVN distribution (for the prior of phi)
estimate_params = function(phi_list) {
  phi_samps = do.call(rbind, lapply(phi_list, getLowerTri))
  mu_phi = apply(phi_samps, 2, mean)
  # sigma_phi = cov(phi_samps)
  sigma_phi = apply(phi_samps, 2, var)
  return(list(mu = mu_phi, sigma = sigma_phi))
}



## log prior functions

lp_phi = function(x, z) {
  p = z$p
  sigma_phi = z$sigma_phi
  mu_phi = z$mu_phi
  c(-p/2 * log(2*pi) - 1/2 * log(det(sigma_phi)) -
      1/2 * t(x-mu_phi) %*% solve(sigma_phi) %*% (x-mu_phi))
}


u = phi_samps[1,]

sum(dnorm(u, mean = z0$mu_phi, sd = sqrt(diag(z0$sigma_phi)), log = TRUE))
lp_phi(u, z0)

lp_phi = function(x, z) {
  p = z$p
  sigma_phi = z$sigma_phi
  mu_phi = z$mu_phi
  tmp = 0
  for (i in 1:p) {
    tmp = tmp - 1/2 * log(2 * pi * sigma_phi[i]) -
      1/(2 * sigma_phi[i]) * (x[i] - mu_phi[i])^2
  }
  return(tmp)
}


logpost_phi = function(x, z) {
  loglik(x, z) + lp_phi(x, z)
}
grad_logpost_phi = function(x, z) {
  # grad_loglik(x, z) + grad_lp_phi(x, z)
  grad_loglik(x, z) + pracma::grad(lp_phi, x0 = x, z = z)
}
hess_logpost_phi = function(x, z) {
  # hess_loglik(x, z) + hess_lp_phi(x, z)
  hess_loglik(x, z) + pracma::hessian(lp_phi, x0 = x, z = z)
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


