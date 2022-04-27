

getLowerTri = function(x) {
  x[lower.tri(x, diag = TRUE)]
}


## sample_beta(): generate PRIOR samples of beta
sample_beta = function(m, k) {
  beta = matrix(0, m, k)
  for (i in 1:m) {
    for (j in 1:min(k, i)) {
      if (i == j) {
        # sample from chi distribution
        beta[i,i] = chi::rchi(1, df = k-i+1)
      } else {
        # sample from normal distribution
        beta[i,j] = rnorm(1, 0, C0)
      }
    }
  }
  return(beta)
}


## sample_sigma(): generate PRIOR samples of sigma (via omega, beta)
sample_sigma = function(J, m, k) {
  omega_samps = matrix(
    pscl::rigamma(m * J, alpha = nu / 2, beta = nu * s^2 / 2), J, m
  )

  sigma_list = vector(mode = "list", length = J)
  for (j in 1:J) {
    omega = diag(omega_samps[j,])
    beta = sample_beta(m, k)
    sigma = omega + beta %*% t(beta)
    sigma_list[[j]] = solve(sigma) ## storing sigma^{-1} b/c we take cholesky
  }
  return(sigma_list)
}


## getLowerChol(): compute phi, the LOWER cholesky factor of sigma = phi phi'
getLowerChol = function(sigma) {
  t(chol(sigma))
}


## using the the prior samples of phi, estimate the mean and covariance of the
## approximating MVN distribution (for the prior of phi)
estimate_params = function(phi_list) {
  phi_samps = do.call(rbind, lapply(phi_list, getLowerTri))
  mu_phi = apply(phi_samps, 2, mean)
  sigma_phi = cov(phi_samps)
  return(list(mu = mu_phi, sigma = sigma_phi))
}


#### functions to extract posterior samples

phi_post = function(obj, J, p) {
  samps = matrix(0, J, p) # store vectorized form of phi in matrix row-wise
  for (j in 1:J) {
    omega = diag(obj$Sigma[[j]])
    beta = obj$beta[[j]]
    sigma = omega + beta %*% t(beta)

    ## parametrization 1: sigma^-1 = phi phi'
    sigma_inv = solve(sigma)
    phi = getLowerChol(sigma_inv)

    samps[j,] = phi[lower.tri(phi, diag = TRUE)]
  }
  return(samps)
}


####### helper functions -------------------------------------------------------

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








